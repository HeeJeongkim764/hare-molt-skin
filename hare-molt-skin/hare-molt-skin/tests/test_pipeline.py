#!/usr/bin/env python3
"""실제 데이터 없이 전체 흐름을 검증하는 테스트.

  python tests/test_pipeline.py

1. 합성 Dryad 자료 생성 → prepare_data.py → run_all.py
2. 심어 둔 효과가 검출되는지 확인 (진피유두·모발·분열 점수는 단계 효과 있음, 멜라닌세포·표피는 없음,
   ASIP 단계 효과는 진피유두 점수로 보정하면 대부분 사라짐)
3. 샘플 이름 네 가지 형식(실제 Dryad 의 A_Brown.rsem.out, A-W, LamAW, hare_A_white)을 모두 읽는지,
   읽을 수 없는 이름이면 템플릿을 만들고 멈추는지
4. 같은 입력으로 두 번 돌리면 verify.py 가 '재현 성공'을 내는지
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable


def run(*args, ok=True):
    r = subprocess.run([PY, *map(str, args)], capture_output=True, text=True, cwd=ROOT)
    if ok and r.returncode != 0:
        raise AssertionError(f"실패: {' '.join(map(str, args))}\n{r.stdout}\n{r.stderr}")
    return r


def prepare(tmp, style, samples=None):
    raw, proc = tmp / f"raw_{style}", tmp / f"proc_{style}"
    run(ROOT / "tests" / "make_synthetic.py", "--out", raw, "--name-style", style, "--n-background", "1500")
    args = [ROOT / "scripts" / "prepare_data.py", "--tar", raw / "Lep_ame_RSEM.out.genes.results.tar",
            "--fasta", raw / "Lep_ame_filtered_skin_transcriptome.fasta",
            "--biomart", raw / "biomart_rabbit.tsv", "--out", proc]
    if samples:
        args += ["--samples", samples]
    return raw, proc, args


def main():
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        print("전체 흐름 테스트")

        # 이름 형식별 샘플 정보 읽기
        for style in ["dryad", "dash", "compact", "words"]:
            raw, proc, args = prepare(tmp, style)
            run(*args)
            got = pd.read_csv(proc / "samples.tsv", sep="\t").sort_values("sample").reset_index(drop=True)
            want = pd.read_csv(raw / "truth_samples.tsv", sep="\t").sort_values("sample").reset_index(drop=True)
            assert got.equals(want), f"{style}: 샘플 정보 추정 오류\n{got}\n{want}"
        print("  샘플 이름 네 형식(실제 Dryad 형식 포함) 자동 인식: 통과")

        raw, proc, args = prepare(tmp, "opaque")
        r = run(*args, ok=False)
        assert r.returncode != 0 and (proc / "samples_template.tsv").exists(), "읽을 수 없는 이름에서 템플릿을 만들지 않음"
        run(*prepare(tmp, "opaque", samples=raw / "truth_samples.tsv")[2])
        print("  이름을 못 읽으면 템플릿 생성 후 --samples 로 진행: 통과")

        # 분석
        proc = tmp / "proc_dash"
        res = tmp / "results"
        run(ROOT / "run_all.py", "--data-dir", proc, "--out-dir", res)
        st = pd.read_csv(res / "tables" / "stage_tests.tsv", sep="\t").set_index("feature")
        for s in ["dermal_papilla", "hair_shaft", "proliferation", "ASIP"]:
            assert st.loc[s, "p_stage"] < 0.01, f"{s}: 심어 둔 단계 효과를 못 찾음 (p={st.loc[s, 'p_stage']})"
        for s in ["melanocyte", "epidermis", "inner_root_sheath"]:
            assert st.loc[s, "p_stage"] > 0.05, f"{s}: 없는 단계 효과를 찾음 (p={st.loc[s, 'p_stage']})"
        assert st.loc["dermal_papilla", "diff_W_minus_B"] > 0 and st.loc["proliferation", "diff_W_minus_B"] < 0
        assert st["n_permutations"].dropna().eq(7776).all(), "순열 수가 7,776 이 아님"
        adj = pd.read_csv(res / "tables" / "asip_adjusted.tsv", sep="\t").set_index("covariate")
        assert adj.loc["dermal_papilla", "fraction_stage_effect_remaining"] < 0.2
        assert adj.loc["melanocyte", "fraction_stage_effect_remaining"] > 0.5
        for fig in ["fig1_sample_correlation", "fig2_pca", "fig3_marker_scores", "fig4_asip_adjusted"]:
            assert (res / "figures" / f"{fig}.png").stat().st_size > 10_000, f"{fig}.png 이 비었음"
        assert "정확 순열검정" in (res / "summary.md").read_text(encoding="utf-8")
        print("  심어 둔 효과 검출, 그림 4개, 요약 생성: 통과")

        # 재현성
        run(ROOT / "scripts" / "verify.py", "--results", res, "--update")
        res2 = tmp / "results2"
        run(ROOT / "run_all.py", "--data-dir", proc, "--out-dir", res2)
        (res2 / "reference").mkdir()
        for p in (res / "reference").glob("*.tsv"):
            (res2 / "reference" / p.name).write_bytes(p.read_bytes())
        out = run(ROOT / "scripts" / "verify.py", "--results", res2).stdout
        assert "재현 성공" in out, out
        print("  두 번 실행 결과 일치 (verify.py): 통과")
    print("모든 테스트 통과")


if __name__ == "__main__":
    main()
