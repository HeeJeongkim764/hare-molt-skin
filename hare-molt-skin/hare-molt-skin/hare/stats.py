"""PCA, 표지 점수, 개체 내 정확 순열검정."""
from itertools import permutations, product

import numpy as np
import pandas as pd

from . import STAGES

TOP_VAR_GENES = 2000
MIN_MARKERS = 2


# ---------------------------------------------------------------- PCA
def pca(logexpr, n_top=TOP_VAR_GENES):
    """변동이 큰 유전자 n_top 개로 PCA (SVD). 반환: 좌표 DataFrame, 설명 분산 비율."""
    var = logexpr.var(axis=1).sort_values(ascending=False)
    x = logexpr.loc[var.index[:n_top]].T.values          # 샘플 × 유전자
    x = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    coords = pd.DataFrame(u[:, :3] * s[:3], index=logexpr.columns, columns=["PC1", "PC2", "PC3"])
    ratio = (s ** 2) / (s ** 2).sum()
    return coords, ratio[:3]


# ---------------------------------------------------------------- 표지 점수
def marker_scores(logexpr, markers):
    """세트마다 유전자별 z-점수(샘플 간)의 평균. 반환: 점수표(세트 × 샘플), 유전자 포함 현황."""
    z = logexpr.sub(logexpr.mean(axis=1), axis=0).div(logexpr.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    scores, coverage = {}, []
    for name, sub in markers.groupby("set", sort=False):
        genes = list(sub["gene"].str.upper())
        found = [g for g in genes if g in z.index and z.loc[g].notna().all()]
        coverage.append({"set": name, "n_defined": len(genes), "n_found": len(found),
                         "genes_found": ",".join(found),
                         "genes_missing": ",".join(g for g in genes if g not in found)})
        scores[name] = z.loc[found].mean(axis=0) if len(found) >= MIN_MARKERS else pd.Series(np.nan, index=z.columns)
    return pd.DataFrame(scores).T, pd.DataFrame(coverage)


# ---------------------------------------------------------------- 정확 순열검정
def to_matrix(values, samples):
    """샘플별 값(Series) → 개체 × 단계(W, I, B) 행렬."""
    m = samples.assign(v=samples["sample"].map(values)).pivot(index="individual", columns="stage", values="v")
    return m[STAGES].values.astype(float)


def _all_relabelings(n_ind):
    """개체마다 단계 라벨 3개를 섞는 모든 경우: (6^n_ind, n_ind, 3) 인덱스 배열."""
    perms = np.array(list(permutations(range(3))))        # (6, 3)
    combos = np.array(list(product(range(6), repeat=n_ind)))
    return perms[combos]


def stage_stats(mat):
    """mat: (..., n_ind, 3) 개체 중심화된 값. 반환: (단계 간 제곱합, W−B 평균 차이)."""
    means = mat.mean(axis=-2)                              # (..., 3)
    ss = mat.shape[-2] * (means ** 2).sum(axis=-1)
    return ss, means[..., 0] - means[..., 2]


def exact_test(mat):
    """개체 안에서만 라벨을 섞는 정확 순열검정 (5개체면 7,776가지 전부).
    mat 은 개체 × 단계 행렬(중심화 전이어도 됨, 내부에서 개체 중심화)."""
    c = mat - mat.mean(axis=1, keepdims=True)
    n_ind = c.shape[0]
    idx = _all_relabelings(n_ind)                          # (P, n_ind, 3)
    permuted = np.take_along_axis(np.broadcast_to(c, idx.shape), idx, axis=2)
    ss_obs, d_obs = stage_stats(c)
    ss_p, d_p = stage_stats(permuted)
    eps = 1e-12
    return {"ss_between": float(ss_obs), "p_stage": float(np.mean(ss_p >= ss_obs - eps)),
            "diff_W_minus_B": float(d_obs), "p_W_vs_B": float(np.mean(np.abs(d_p) >= abs(d_obs) - eps)),
            "n_permutations": int(len(ss_p))}


def residualize(y, x):
    """개체 중심화된 y 를 개체 중심화된 x 로 회귀(절편 없음)한 잔차와 기울기."""
    yc = y - y.mean(axis=1, keepdims=True)
    xc = x - x.mean(axis=1, keepdims=True)
    denom = (xc ** 2).sum()
    beta = (yc * xc).sum() / denom if denom > 0 else 0.0
    return yc - beta * xc, float(beta)


def within_corr(y, x):
    yc = (y - y.mean(axis=1, keepdims=True)).ravel()
    xc = (x - x.mean(axis=1, keepdims=True)).ravel()
    return float(np.corrcoef(yc, xc)[0, 1])
