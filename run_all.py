#!/usr/bin/env python3
"""전체 분석 한 번에 실행.

  python run_all.py

입력 : data/processed/tpm.tsv, samples.tsv (레포에 포함), config/marker_sets.tsv
출력 : results/figures/*.png, results/tables/*.tsv, results/summary.md
무작위 요소가 없어서(정확 순열검정) 어느 컴퓨터에서 돌려도 결과가 같습니다.
"""
import argparse
import platform
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from hare import STAGES, TARGET_GENE, data, figures, stats

ROOT = Path(__file__).resolve().parent


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=ROOT / "data" / "processed")
    ap.add_argument("--markers", default=ROOT / "config" / "marker_sets.tsv")
    ap.add_argument("--out-dir", default=ROOT / "results")
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out_dir)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    # 1. 데이터
    print("[1/5] 데이터 읽기·정규화")
    tpm, samples = data.load(a.data_dir)
    logexpr, n_kept = data.filter_and_log(tpm)
    inds = sorted(samples["individual"].unique())
    print(f"      샘플 {len(samples)}개 (개체 {len(inds)} × 단계 3), 유전자 {len(tpm)}개 중 {n_kept}개 사용")

    # 2. 품질 확인
    print("[2/5] 품질 확인: 샘플 상관, PCA")
    figures.sample_correlation(logexpr, samples, out / "figures" / "fig1_sample_correlation.png")
    raw_pc, raw_ratio = stats.pca(logexpr)
    cen = data.center_by_individual(logexpr, samples)
    cen_pc, cen_ratio = stats.pca(cen)
    figures.pca_panels(raw_pc, raw_ratio, cen_pc, cen_ratio, samples, out / "figures" / "fig2_pca.png")
    pcs = raw_pc.join(cen_pc, rsuffix="_centered").join(samples.set_index("sample")[["individual", "stage"]])
    pcs.round(6).to_csv(out / "tables" / "pca_coordinates.tsv", sep="\t", index_label="sample")
    corr = np.corrcoef(logexpr.T.values)
    same_ind = np.equal.outer(samples["individual"].values, samples["individual"].values)
    same_stage = np.equal.outer(samples["stage"].values, samples["stage"].values)
    off = ~np.eye(len(samples), dtype=bool)
    qc = {"mean_r_same_individual": corr[same_ind & off].mean(),
          "mean_r_same_stage_diff_individual": corr[same_stage & ~same_ind].mean(),
          "mean_r_diff_both": corr[~same_stage & ~same_ind].mean(),
          "pc1_var_raw": raw_ratio[0], "pc1_var_centered": cen_ratio[0]}

    # 3. 세포 유형 표지 점수
    print("[3/5] 세포 유형 표지 점수")
    markers = pd.read_csv(a.markers, sep="\t", dtype=str)
    scores, coverage = stats.marker_scores(logexpr, markers)
    coverage.to_csv(out / "tables" / "marker_coverage.tsv", sep="\t", index=False)
    scores.T.join(samples.set_index("sample")[["individual", "stage"]]).round(6) \
        .to_csv(out / "tables" / "marker_scores.tsv", sep="\t", index_label="sample")
    score_mats, tests, rows = {}, {}, []
    for name in scores.index:
        mat = stats.to_matrix(scores.loc[name], samples)
        score_mats[name] = mat
        if np.isnan(mat).any():
            rows.append({"feature": name, "type": "marker_set", "note": "표지 유전자 부족"})
            continue
        t = stats.exact_test(mat)
        tests[name] = t
        rows.append({"feature": name, "type": "marker_set", **t})
    labels = {r.set: f"{r.set.replace('_', ' ')} ({r.n_found}/{r.n_defined} genes)" for r in coverage.itertuples()}
    figures.marker_score_panels(score_mats, tests, inds, samples, out / "figures" / "fig3_marker_scores.png", labels)

    # 4. 어구티
    print(f"[4/5] {TARGET_GENE} 단계 효과와 세포 구성 보정")
    if TARGET_GENE not in logexpr.index:
        sys.exit(f"{TARGET_GENE} 가 필터 후 발현표에 없습니다 (발현이 낮거나 주석 누락). data/processed/tpm.tsv 확인")
    asip = stats.to_matrix(logexpr.loc[TARGET_GENE], samples)
    asip_test = stats.exact_test(asip)
    rows.append({"feature": TARGET_GENE, "type": "gene", **asip_test})
    adj = []
    for name, t in tests.items():
        resid, beta = stats.residualize(asip, score_mats[name])
        at = stats.exact_test(resid)
        adj.append({"covariate": name, "within_individual_r": stats.within_corr(asip, score_mats[name]),
                    "beta": beta, "ss_before": asip_test["ss_between"], "ss_after": at["ss_between"],
                    "fraction_stage_effect_remaining": at["ss_between"] / asip_test["ss_between"]
                    if asip_test["ss_between"] > 0 else np.nan,
                    "p_stage_after": at["p_stage"], "diff_W_minus_B_after": at["diff_W_minus_B"],
                    "p_W_vs_B_after": at["p_W_vs_B"]})
    pd.DataFrame(rows).round(6).to_csv(out / "tables" / "stage_tests.tsv", sep="\t", index=False)
    adj_df = pd.DataFrame(adj)
    adj_df.round(6).to_csv(out / "tables" / "asip_adjusted.tsv", sep="\t", index=False)
    figures.asip_panels(asip, inds, samples, adj, asip_test, out / "figures" / "fig4_asip_adjusted.png")

    # 5. 요약
    print("[5/5] 요약 작성")
    elapsed = time.time() - t0
    write_summary(out, samples, tpm, n_kept, qc, coverage, rows, asip_test, adj_df, elapsed)
    print(f"\n완료 ({elapsed:.1f}초) → {out / 'summary.md'}")


def write_summary(out, samples, tpm, n_kept, qc, coverage, rows, asip_test, adj_df, elapsed):
    f = lambda x: f"{x:.3f}"
    L = ["# 분석 결과 요약 (run_all.py 가 자동 생성)", "",
         "## 무엇을 했나", "",
         f"- 샘플 {len(samples)}개: 개체 {samples['individual'].nunique()}마리 × 털갈이 단계 3개(W 흰색, I 중간, B 갈색)",
         f"- 유전자 {len(tpm)}개 중 TPM ≥ {data.MIN_TPM} 인 샘플이 {data.MIN_SAMPLES}개 이상인 {n_kept}개 사용, log2(TPM+1)",
         "- 품질 확인: 샘플 간 상관, PCA (원래 값 / 개체 평균 제거 후)",
         "- 세포 유형 표지 점수: 세트별 유전자 z-점수 평균",
         f"- {TARGET_GENE} 단계 효과를 각 표지 점수로 보정했을 때 얼마나 남는지",
         f"- 통계: 개체 안에서만 단계 라벨을 섞는 정확 순열검정 ({asip_test['n_permutations']}가지 전부)", "",
         "## 품질 확인", "",
         f"- 같은 개체 샘플끼리 평균 상관 r = {f(qc['mean_r_same_individual'])}",
         f"- 같은 단계·다른 개체끼리 r = {f(qc['mean_r_same_stage_diff_individual'])}",
         f"- 단계·개체 모두 다를 때 r = {f(qc['mean_r_diff_both'])}",
         f"- PC1 설명 분산: 원래 {qc['pc1_var_raw'] * 100:.1f}% → 개체 평균 제거 후 {qc['pc1_var_centered'] * 100:.1f}%", "",
         "## 표지 유전자 포함 현황", "", "| 세트 | 정의 | 찾음 | 없는 유전자 |", "|---|---|---|---|"]
    for r in coverage.itertuples():
        L.append(f"| {r.set} | {r.n_defined} | {r.n_found} | {r.genes_missing or '-'} |")
    weak = [r.set for r in coverage.itertuples() if r.n_found < 0.5 * r.n_defined]
    if weak:
        L += ["", f"> 주의: {', '.join(weak)} 세트는 정의한 유전자의 절반 미만만 찾아서 점수의 신뢰도가 낮습니다."]
    L += ["", "## 단계 효과 (정확 순열검정)", "", "| 특징 | 단계 간 제곱합 | p (단계) | W−B 차이 | p (W vs B) |",
          "|---|---|---|---|---|"]
    for r in rows:
        if "p_stage" in r:
            L.append(f"| {r['feature']} | {f(r['ss_between'])} | {r['p_stage']:.4f} | {f(r['diff_W_minus_B'])} | {r['p_W_vs_B']:.4f} |")
        else:
            L.append(f"| {r['feature']} | - | - | - | {r.get('note', '')} |")
    L += ["", f"## {TARGET_GENE} 를 세포 구성 점수로 보정", "",
          "| 보정 변수 | 개체 내 상관 r | 남은 단계 효과 비율 | p (보정 후) |", "|---|---|---|---|"]
    for r in adj_df.itertuples():
        L.append(f"| {r.covariate} | {f(r.within_individual_r)} | {f(r.fraction_stage_effect_remaining)} | {r.p_stage_after:.4f} |")
    min_p = 6 / asip_test["n_permutations"]
    L += ["", "## 해석할 때 주의", "",
          f"- 개체가 5마리라 가능한 가장 작은 p 는 {min_p:.5f} 입니다 (모든 개체에서 라벨을 똑같이 바꾸는 6가지는 같은 통계량).",
          "- 표지 점수는 세포 비율을 직접 잰 값이 아니라 표지 유전자 발현의 평균입니다.",
          "- 보정 변수도 단계에 따라 변하므로, 보정 후 효과가 줄었다는 것은 \"구성 변화와 구분이 안 된다\"는 뜻이지 원인이 증명된 것은 아닙니다.",
          "- 단계를 따라 움직이는 점수라면 무엇으로 보정해도 효과가 줄 수 있습니다. 어느 세포 유형 점수가 ASIP 와 특히 잘 맞는지(개체 내 상관)를 함께 보고 해석하세요.", "",
          "## 실행 환경", "",
          f"- Python {platform.python_version()}, pandas {pd.__version__}, numpy {np.__version__}, matplotlib {matplotlib.__version__}",
          f"- 실행 시간 {elapsed:.1f}초", ""]
    (out / "summary.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
