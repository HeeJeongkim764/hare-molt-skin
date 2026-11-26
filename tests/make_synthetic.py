#!/usr/bin/env python3
"""Dryad 원자료와 같은 형식의 합성 데이터 생성 (실제 데이터 없이 전체 흐름 검증용).

만드는 것: RSEM *.genes.results 15개를 담은 tar, Trinity FASTA(헤더에 Ensembl ID), BioMart 표.
심어 둔 효과 (개체 중심화 후):
  - dermal_papilla 점수: W 에서 높고 B 에서 낮음
  - ASIP: dermal_papilla 와 같이 움직임 (W 높음) → 진피유두 점수로 보정하면 단계 효과가 크게 줄어야 함
  - proliferation, hair_shaft: B 에서 높음 (anagen)
  - epidermis: 단계 효과 없음
  - 개체 간 차이를 크게 넣어, 개체 중심화가 필요한 상황을 재현
"""
import argparse
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STAGE_EFFECT = {"W": 1.0, "I": 0.0, "B": -1.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-background", type=int, default=3000)
    ap.add_argument("--name-style", default="dash", choices=["dash", "compact", "words", "opaque", "dryad"])
    a = ap.parse_args()
    rng = np.random.default_rng(7)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    markers = pd.read_csv(ROOT / "config" / "marker_sets.tsv", sep="\t")
    genes = list(dict.fromkeys(list(markers["gene"]) + ["ASIP"] + [f"GENE{i:05d}" for i in range(a.n_background)]))
    gset = dict(zip(markers["gene"], markers["set"]))
    n = len(genes)
    base = rng.uniform(1, 8, n)                                   # log2 기준 발현
    ind_eff = {i: rng.normal(0, 0.8, n) for i in "ABCDE"}          # 큰 개체 효과

    samples, cols = [], {}
    for ind in "ABCDE":
        for st in "WIB":
            e = STAGE_EFFECT[st]
            lg = base + ind_eff[ind] + rng.normal(0, 0.15, n)
            for k, g in enumerate(genes):
                s = gset.get(g)
                if s == "dermal_papilla":
                    lg[k] += 0.9 * e
                elif s in ("proliferation", "hair_shaft"):
                    lg[k] -= 0.8 * e
                elif g == "ASIP":
                    lg[k] += 1.5 * e
            name = {"dash": f"{ind}-{st}", "compact": f"Lam{ind}{st}", "words": f"hare_{ind}_{ {'W': 'white', 'I': 'intermediate', 'B': 'brown'}[st] }",
                    "opaque": f"S{len(samples) + 1:02d}",
                    "dryad": f"{ind}_{ {'W': 'White', 'I': 'Intermediate', 'B': 'Brown'}[st] }"}[a.name_style]
            cols[name] = np.maximum(2 ** lg - 1, 0)
            samples.append((name, ind, st))

    # Trinity 조각: 유전자 하나에 조각 1~2개 (합산 검증)
    comps, comp_gene = [], []
    for k, g in enumerate(genes):
        for c in range(1 if k % 7 else 2):
            comps.append(f"comp{k}_c{c}")
            comp_gene.append(k)
    ens = {k: f"ENSOCUG{k:011d}" for k in range(n)}

    tar_path = out / "Lep_ame_RSEM.out.genes.results.tar"
    with tarfile.open(tar_path, "w") as tf:
        for name, vals in cols.items():
            frac = np.array([0.6 if (k % 7 == 0) else 1.0 for k in comp_gene])
            frac = np.where([(k % 7 == 0) and c.endswith("c1") for k, c in zip(comp_gene, comps)], 0.4, frac)
            tpm = vals[comp_gene] * frac
            df = pd.DataFrame({"gene_id": comps, "transcript_id(s)": [c + "_seq1" for c in comps],
                               "length": 1500, "effective_length": 1300,
                               "expected_count": np.round(tpm * 3, 2), "TPM": np.round(tpm, 4),
                               "FPKM": np.round(tpm * 0.9, 4)})
            b = df.to_csv(sep="\t", index=False).encode()
            suffix = ".rsem.out.genes.results" if a.name_style == "dryad" else ".genes.results"
            ti = tarfile.TarInfo(f"{name}{suffix}" if a.name_style == "dryad" else f"RSEM_results/{name}{suffix}")
            ti.size = len(b)
            tf.addfile(ti, io.BytesIO(b))

    with open(out / "Lep_ame_filtered_skin_transcriptome.fasta", "w") as f:
        for c, k in zip(comps, comp_gene):
            f.write(f">{c}_seq1 len=1500 {ens[k]} best_hit\nACGTACGTACGT\n")
    pd.DataFrame({"Gene stable ID": [ens[k] for k in range(n)], "Gene name": genes}) \
        .to_csv(out / "biomart_rabbit.tsv", sep="\t", index=False)
    pd.DataFrame(samples, columns=["sample", "individual", "stage"]).to_csv(out / "truth_samples.tsv", sep="\t", index=False)
    print(f"합성 데이터 → {out}")


if __name__ == "__main__":
    main()
