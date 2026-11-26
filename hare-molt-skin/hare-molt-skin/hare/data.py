"""데이터 읽기, 필터, 정규화."""
from pathlib import Path

import numpy as np
import pandas as pd

from . import STAGES

MIN_TPM = 1.0       # 발현으로 볼 최소 TPM
MIN_SAMPLES = 3     # 그 이상인 샘플 수 (한 단계 = 5개 중 최소 3개 정도)


def load(data_dir):
    d = Path(data_dir)
    if not (d / "tpm.tsv").exists() or not (d / "samples.tsv").exists():
        raise SystemExit(f"{d} 에 tpm.tsv / samples.tsv 가 없습니다.\n"
                         "처음 한 번은 README 의 '데이터 준비' 대로 scripts/prepare_data.py 를 실행하세요.")
    tpm = pd.read_csv(d / "tpm.tsv", sep="\t", index_col=0)
    samples = pd.read_csv(d / "samples.tsv", sep="\t", dtype=str)
    tpm.index = tpm.index.astype(str).str.upper()
    samples["stage"] = samples["stage"].str.upper()
    check_design(samples)
    # 개체 → 단계(W, I, B) 순서로 정렬
    samples["stage_order"] = samples["stage"].map({s: i for i, s in enumerate(STAGES)})
    samples = samples.sort_values(["individual", "stage_order"]).drop(columns="stage_order").reset_index(drop=True)
    missing = set(samples["sample"]) - set(tpm.columns)
    if missing:
        raise SystemExit(f"tpm.tsv 에 없는 샘플: {sorted(missing)}")
    return tpm[samples["sample"]], samples


def check_design(samples):
    bad_stage = set(samples["stage"]) - set(STAGES)
    if bad_stage:
        raise SystemExit(f"stage 는 W/I/B 여야 합니다: {bad_stage}")
    tab = samples.groupby(["individual", "stage"]).size().unstack(fill_value=0)
    if not ((tab.reindex(columns=STAGES, fill_value=0) == 1).all().all()):
        raise SystemExit(f"개체마다 W/I/B 가 정확히 하나씩 있어야 합니다:\n{tab}")


def filter_and_log(tpm):
    keep = (tpm >= MIN_TPM).sum(axis=1) >= MIN_SAMPLES
    return np.log2(tpm.loc[keep] + 1), int(keep.sum())


def center_by_individual(values, samples):
    """개체별 평균을 빼서 개체 차이를 제거 (행: 특징, 열: 샘플)."""
    out = values.copy()
    for ind, sub in samples.groupby("individual"):
        cols = list(sub["sample"])
        out[cols] = values[cols].sub(values[cols].mean(axis=1), axis=0)
    return out
