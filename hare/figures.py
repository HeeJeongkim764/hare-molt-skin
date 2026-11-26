"""그림 (정적 PNG). 어느 컴퓨터에서도 글꼴 문제가 없도록 그림 안 글자는 영어로."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import STAGES, STAGE_NAMES, TARGET_GENE  # noqa: E402

STAGE_COLORS = {"W": "#257de0", "I": "#00b894", "B": "#f47230"}   # 검증된 범주형 팔레트 앞 3칸
IND_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
TEXT, TEXT2, GRID = "#191919", "#666666", "#e6e6e6"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": "#bdbdbd", "axes.labelcolor": TEXT, "xtick.color": TEXT2,
    "ytick.color": TEXT2, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
})


def _ind_marker(samples):
    inds = sorted(samples["individual"].unique())
    return {ind: IND_MARKERS[i % len(IND_MARKERS)] for i, ind in enumerate(inds)}


def sample_correlation(logexpr, samples, path):
    corr = np.corrcoef(logexpr.T.values)
    labels = [f"{r.individual}-{r.stage}" for r in samples.itertuples()]
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    im = ax.imshow(corr, cmap="Blues", vmin=np.floor(corr.min() * 100) / 100, vmax=1)
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    for k in range(3, len(labels), 3):                    # 개체 경계
        ax.axhline(k - 0.5, color="white", lw=2)
        ax.axvline(k - 0.5, color="white", lw=2)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Pearson r (log2 TPM+1)")
    ax.set_title("Sample-to-sample correlation (grouped by individual)")
    ax.spines[:].set_visible(False)
    fig.savefig(path)
    plt.close(fig)


def pca_panels(raw, raw_ratio, cen, cen_ratio, samples, path):
    mk = _ind_marker(samples)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, coords, ratio, title in [(axes[0], raw, raw_ratio, "PCA: log2 expression"),
                                     (axes[1], cen, cen_ratio, "PCA: individual means removed")]:
        for r in samples.itertuples():
            ax.scatter(coords.loc[r.sample, "PC1"], coords.loc[r.sample, "PC2"], s=70,
                       color=STAGE_COLORS[r.stage], marker=mk[r.individual], edgecolor="white", linewidth=1.5, zorder=3)
        ax.set_xlabel(f"PC1 ({ratio[0] * 100:.1f}%)")
        ax.set_ylabel(f"PC2 ({ratio[1] * 100:.1f}%)")
        ax.set_title(title)
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    stage_h = [plt.Line2D([], [], marker="o", ls="", color=STAGE_COLORS[s], ms=8, label=STAGE_NAMES[s]) for s in STAGES]
    ind_h = [plt.Line2D([], [], marker=m, ls="", color=TEXT2, ms=7, label=f"hare {i}") for i, m in mk.items()]
    fig.legend(handles=stage_h + ind_h, loc="lower center", ncol=len(stage_h) + len(ind_h), bbox_to_anchor=(0.5, -0.13))
    fig.savefig(path)
    plt.close(fig)


def _paired(ax, mat, inds, mk, ylabel):
    x = np.arange(3)
    for row, ind in zip(mat, inds):
        ax.plot(x, row, color="#b0b0b0", lw=1.2, zorder=1)
        for j, s in enumerate(STAGES):
            ax.scatter(j, row[j], color=STAGE_COLORS[s], marker=mk[ind], s=36, edgecolor="white", linewidth=1, zorder=3)
    ax.plot(x, np.nanmean(mat, axis=0), color=TEXT, lw=2, zorder=2)
    ax.set_xticks(x, [STAGE_NAMES[s] for s in STAGES])
    ax.set_xlim(-0.35, 2.35)
    ax.axhline(0, color=GRID, lw=0.8, zorder=0)
    ax.set_ylabel(ylabel)


def marker_score_panels(score_mats, tests, inds, samples, path, labels=None):
    mk = _ind_marker(samples)
    names = list(score_mats)
    ncol = 3
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 3.1 * nrow), squeeze=False)
    for ax, name in zip(axes.flat, names):
        _paired(ax, score_mats[name], inds, mk, "score (mean z)")
        t = tests.get(name)
        sub = f"stage p = {t['p_stage']:.4f}" if t else "too few marker genes"
        label = (labels or {}).get(name, name.replace("_", " "))
        ax.set_title(f"{label}\n{sub}", fontsize=10)
    for ax in list(axes.flat)[len(names):]:
        ax.axis("off")
    fig.suptitle("Cell-type marker scores across the spring molt (grey line = one hare, black = mean)", y=1.0)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def asip_panels(asip_mat, inds, samples, adj, raw_test, path):
    mk = _ind_marker(samples)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 1.3]})
    centered = asip_mat - asip_mat.mean(axis=1, keepdims=True)
    _paired(axes[0], centered, inds, mk, f"{TARGET_GENE} log2 expression\n(individual mean removed)")
    axes[0].set_title(f"{TARGET_GENE} across the molt\nstage p = {raw_test['p_stage']:.4f}")
    ax = axes[1]
    labels = ["none"] + [a["covariate"].replace("_", " ") for a in adj]
    vals = [1.0] + [a["fraction_stage_effect_remaining"] for a in adj]
    pv = [raw_test["p_stage"]] + [a["p_stage_after"] for a in adj]
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, height=0.55, color=["#4d4d4d"] + ["#257de0"] * len(adj))
    for yi, v, p in zip(y, vals, pv):
        ax.text(v + 0.02, yi, f"{v:.2f}  (p = {p:.4f})", va="center", fontsize=9, color=TEXT)
    ax.set_yticks(y, labels)
    ax.set_xlim(0, max(1.45, max(vals) + 0.5))
    ax.set_xlabel(f"{TARGET_GENE} stage effect remaining after adjustment (1 = unchanged)")
    ax.set_title(f"Does cell composition explain the {TARGET_GENE} change?")
    ax.axvline(1, color=GRID, lw=0.8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
