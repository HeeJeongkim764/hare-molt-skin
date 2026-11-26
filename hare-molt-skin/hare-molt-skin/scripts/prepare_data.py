#!/usr/bin/env python3
"""Dryad 원자료 → 분석용 발현표 (처음 한 번만 실행).

입력 (Dryad doi:10.5061/dryad.7v1m1 에서 브라우저로 내려받아 data/raw/ 에 두기)
  필수  Lep_ame_RSEM.out.genes.results.tar   15개 샘플의 RSEM 유전자 단위 결과
  선택  Lep_ame_filtered_skin_transcriptome.fasta   Trinity 조각 → Ensembl ID 주석이 헤더에 있을 때
  선택  Ensembl BioMart 표 (Gene stable ID, [Transcript stable ID,] Gene name)   ID → 유전자 기호

출력 (data/processed/, 레포에 커밋)
  tpm.tsv       유전자 기호 × 샘플 TPM (같은 기호의 조각은 합산)
  counts.tsv    유전자 기호 × 샘플 expected_count
  samples.tsv   sample, individual(A–E), stage(W/I/B)
  prepare_report.txt   무엇을 어떻게 처리했는지 기록

사용 예
  python scripts/prepare_data.py --tar data/raw/Lep_ame_RSEM.out.genes.results.tar \
      --fasta data/raw/Lep_ame_filtered_skin_transcriptome.fasta \
      --biomart data/raw/biomart_rabbit.tsv data/raw/biomart_mouse.tsv
샘플 이름에서 개체/단계를 자동으로 못 읽으면 data/processed/samples_template.tsv 를 만들고 멈춥니다.
그 파일을 채운 뒤 --samples 로 넘기세요.
"""
import argparse
import csv
import io
import re
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"

RE_ENS_GENE = re.compile(r"\bENS[A-Z]{0,6}G\d{11}\b")
RE_ENS_TX = re.compile(r"\bENS[A-Z]{0,6}T\d{11}\b")
RE_ENS_PROT = re.compile(r"\bENS[A-Z]{0,6}P\d{11}\b")
RE_TRINITY_SUFFIX = re.compile(r"_(seq|i)\d+$")
STAGE_WORDS = {"W": "W", "WHITE": "W", "I": "I", "INT": "I", "INTER": "I", "INTERMEDIATE": "I",
               "B": "B", "BROWN": "B"}
INDIVIDUALS = set("ABCDE")


# ---------------------------------------------------------------- RSEM 읽기
def read_rsem_tar(path):
    """tar 안의 *.genes.results 를 {샘플: DataFrame} 으로."""
    out = {}
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            name = Path(m.name).name
            if not m.isfile() or name.startswith(".") or "genes.results" not in name:
                continue
            raw = tf.extractfile(m).read()
            if name.endswith(".gz"):
                import gzip
                raw = gzip.decompress(raw)
            df = pd.read_csv(io.BytesIO(raw), sep="\t")
            sample = re.sub(r"\.genes\.results(\.gz)?$", "", name)
            sample = re.sub(r"(_|\.)?RSEM(\.out)?$", "", sample, flags=re.I)
            out[sample] = df
    return out


def read_rsem_dir(path):
    out = {}
    for p in sorted(Path(path).glob("*genes.results*")):
        sample = re.sub(r"\.genes\.results(\.gz)?$", "", p.name)
        out[sample] = pd.read_csv(p, sep="\t")
    return out


# ---------------------------------------------------------------- 샘플 정보
def parse_sample(name):
    """샘플 이름에서 (개체, 단계) 추정. 실패하면 (None, None)."""
    up = name.upper()
    # 1) 구분자가 있으면 토큰으로: 단계는 뒤에서부터, 개체는 남은 토큰 중 한 글자 A–E
    toks = [t for t in re.split(r"[^A-Z0-9]+", up) if t]
    if len(toks) >= 2:
        si = next((i for i in range(len(toks) - 1, -1, -1) if toks[i] in STAGE_WORDS), None)
        if si is not None:
            ind = next((t for i, t in enumerate(toks) if i != si and t in INDIVIDUALS), None)
            if ind:
                return ind, STAGE_WORDS[toks[si]]
    # 2) 붙어 있는 이름(예: LamAW, AW1)은 끝 또는 처음의 두 글자만 인정
    compact = re.sub(r"[^A-Z0-9]", "", up)
    m = re.search(r"([A-E])(W|I|B)\d*$", compact) or re.match(r"([A-E])(W|I|B)", compact)
    if m:
        return m.group(1), m.group(2)
    return None, None


def build_samples(names, samples_file, out_dir):
    if samples_file:
        s = pd.read_csv(samples_file, sep="\t", dtype=str)
        need = {"sample", "individual", "stage"}
        if not need <= set(s.columns):
            sys.exit(f"--samples 파일에 {need} 열이 필요합니다")
        missing = set(names) - set(s["sample"])
        if missing:
            sys.exit(f"--samples 에 없는 샘플: {sorted(missing)}")
        return s[s["sample"].isin(names)][["sample", "individual", "stage"]].reset_index(drop=True)
    rows = [(n, *parse_sample(n)) for n in names]
    s = pd.DataFrame(rows, columns=["sample", "individual", "stage"])
    design = s.dropna().groupby(["individual", "stage"]).size()
    ok = (s.notna().all().all() and len(s) == 15 and len(design) == 15 and (design == 1).all())
    if not ok:
        out_dir.mkdir(parents=True, exist_ok=True)
        tmpl = out_dir / "samples_template.tsv"
        s.fillna("").to_csv(tmpl, sep="\t", index=False)
        sys.exit(f"샘플 이름에서 개체(A–E)/단계(W/I/B)를 확실히 읽지 못했습니다.\n"
                 f"{tmpl} 의 individual, stage 열을 채운 뒤 --samples {tmpl} 로 다시 실행하세요.\n"
                 f"(논문 Table S1 또는 SRA 메타데이터 SRR5468532–SRR5468561 참고)")
    return s


# ---------------------------------------------------------------- 유전자 주석
def load_biomart(paths):
    """BioMart 표들 → (Ensembl 유전자 ID → 기호, 전사체·단백질 ID → 유전자 ID)."""
    g2s, t2g = {}, {}
    for p in paths or []:
        df = pd.read_csv(p, sep="\t", dtype=str)
        cols = {c.lower(): c for c in df.columns}
        gcol = cols.get("gene stable id") or cols.get("ensembl_gene_id")
        ncol = cols.get("gene name") or cols.get("external_gene_name")
        tcol = cols.get("transcript stable id") or cols.get("ensembl_transcript_id")
        pcol = cols.get("protein stable id") or cols.get("ensembl_peptide_id")
        if not gcol or not ncol:
            sys.exit(f"{p}: 'Gene stable ID' 와 'Gene name' 열이 필요합니다")
        for _, r in df.iterrows():
            g, n = r[gcol], r[ncol]
            if isinstance(g, str) and isinstance(n, str) and n.strip():
                g2s.setdefault(g.split(".")[0], n.strip().upper())
            for c in (tcol, pcol):
                if c and isinstance(r[c], str) and isinstance(g, str):
                    t2g[r[c].split(".")[0]] = g.split(".")[0]
    return g2s, t2g


def fasta_annotation(path, t2g, g2s):
    """FASTA 헤더에서 Trinity 유전자 ID → Ensembl 유전자 ID.
    헤더에 ID가 여럿이면(예: ENSMUSG…|ENSOCUG…) 논문처럼 토끼(ENSOCUG)를 우선하되,
    유전자 기호가 있는 ID를 고른다 (토끼 ID에 기호가 없으면 마우스 ID 사용)."""
    comp2ens = {}
    with open(path, errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            head = line[1:].strip()
            tx = re.split(r"[\s|]", head, maxsplit=1)[0]
            comp = RE_TRINITY_SUFFIX.sub("", tx)
            if comp in comp2ens:
                continue
            cands = RE_ENS_GENE.findall(head)
            for rx in (RE_ENS_TX, RE_ENS_PROT):          # 전사체·단백질 ID 만 있는 경우
                cands += [t2g[t] for t in rx.findall(head) if t in t2g]
            if not cands:
                continue
            cands.sort(key=lambda g: (0 if g.startswith("ENSOCUG") else 1))
            named = [g for g in cands if g in g2s]
            comp2ens[comp] = named[0] if named else cands[0]
    return comp2ens


def annotate(gene_ids, args, report):
    g2s, t2g = load_biomart(args.biomart)
    if args.annotation:
        a = pd.read_csv(args.annotation, sep="\t", dtype=str)
        m = dict(zip(a.iloc[:, 0], a.iloc[:, 1].str.upper()))
        report.append(f"주석: --annotation 표 사용 ({len(m)}개)")
        return {g: m.get(g) for g in gene_ids}
    n_ens = sum(bool(RE_ENS_GENE.fullmatch(g.split(".")[0])) for g in gene_ids)
    if n_ens > 0.5 * len(gene_ids):
        report.append(f"주석: gene_id 가 Ensembl 유전자 ID ({n_ens}/{len(gene_ids)}) → BioMart 로 기호 변환")
        if not g2s:
            sys.exit("gene_id 가 Ensembl ID 입니다. --biomart 로 Ensembl BioMart 표(Gene stable ID, Gene name)를 주세요")
        return {g: g2s.get(g.split(".")[0]) for g in gene_ids}
    if args.fasta:
        comp2ens = fasta_annotation(args.fasta, t2g, g2s)
        report.append(f"주석: FASTA 헤더에서 Trinity → Ensembl {len(comp2ens)}개 추출")
        if not comp2ens:
            sys.exit("FASTA 헤더에서 Ensembl ID 를 찾지 못했습니다. 헤더 몇 줄을 확인하고 --annotation 표를 직접 만드세요")
        if not g2s:
            sys.exit("Ensembl ID → 유전자 기호 변환에 --biomart 표가 필요합니다")
        return {g: g2s.get(comp2ens.get(g, "")) for g in gene_ids}
    sys.exit("유전자 기호를 붙일 방법이 없습니다. --fasta + --biomart 또는 --annotation 을 주세요.\n"
             f"gene_id 예시: {list(gene_ids)[:5]}")


# ---------------------------------------------------------------- main
def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--tar", help="Lep_ame_RSEM.out.genes.results.tar")
    src.add_argument("--rsem-dir", help="tar 를 이미 푼 폴더")
    ap.add_argument("--fasta", help="Lep_ame_filtered_skin_transcriptome.fasta (선택)")
    ap.add_argument("--biomart", nargs="*", help="Ensembl BioMart TSV (토끼, 마우스)")
    ap.add_argument("--annotation", help="직접 만든 gene_id<TAB>symbol 표 (선택)")
    ap.add_argument("--samples", help="sample, individual, stage 열을 가진 TSV (선택)")
    ap.add_argument("--out", default=str(OUT), help="출력 폴더 (기본 data/processed)")
    args = ap.parse_args()

    report = []
    data = read_rsem_tar(args.tar) if args.tar else read_rsem_dir(args.rsem_dir)
    if not data:
        sys.exit("RSEM *.genes.results 파일을 찾지 못했습니다")
    names = sorted(data)
    report.append(f"RSEM 파일 {len(names)}개: {', '.join(names)}")
    first = data[names[0]]
    for col in ("gene_id", "expected_count", "TPM"):
        if col not in first.columns:
            sys.exit(f"RSEM 파일에 '{col}' 열이 없습니다. 열: {list(first.columns)}")

    out = Path(args.out)
    samples = build_samples(names, args.samples, out)
    gene_ids = list(first["gene_id"].astype(str))
    sym = annotate(gene_ids, args, report)
    n_sym = sum(1 for v in sym.values() if v)
    report.append(f"기호가 붙은 유전자 단위: {n_sym}/{len(gene_ids)}")
    if n_sym < 0.3 * len(gene_ids):
        print(f"[경고] 기호가 붙은 비율이 낮습니다 ({n_sym}/{len(gene_ids)}). 주석 입력을 확인하세요.", file=sys.stderr)

    tpm = pd.DataFrame({s: data[s].set_index("gene_id")["TPM"] for s in samples["sample"]})
    cnt = pd.DataFrame({s: data[s].set_index("gene_id")["expected_count"] for s in samples["sample"]})
    tpm.index = tpm.index.astype(str)
    cnt.index = cnt.index.astype(str)
    keep = [g for g in tpm.index if sym.get(g)]
    symbols = pd.Series({g: sym[g] for g in keep})
    tpm_s = tpm.loc[keep].groupby(symbols).sum().sort_index()
    cnt_s = cnt.loc[keep].groupby(symbols).sum().sort_index()
    report.append(f"같은 기호 합산 후 유전자 수: {len(tpm_s)}")

    out.mkdir(parents=True, exist_ok=True)
    tpm_s.round(4).to_csv(out / "tpm.tsv", sep="\t", index_label="gene")
    cnt_s.round(2).to_csv(out / "counts.tsv", sep="\t", index_label="gene")
    samples.to_csv(out / "samples.tsv", sep="\t", index=False)
    (out / "prepare_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    print(f"\n완료 → {out}")


if __name__ == "__main__":
    main()
