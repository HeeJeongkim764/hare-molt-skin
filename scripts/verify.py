#!/usr/bin/env python3
"""재현성 확인: 방금 만든 results/tables 가 레포에 저장된 기준 결과(results/reference)와 같은지 비교.

  python scripts/verify.py            # 비교
  python scripts/verify.py --update   # 현재 결과를 기준으로 저장 (실제 데이터로 처음 돌린 뒤 한 번)
숫자는 상대 오차 1e-6 까지 같으면 일치로 봅니다 (운영체제·라이브러리 버전에 따른 마지막 자리 차이 허용).
"""
import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RTOL, ATOL = 1e-6, 1e-9


def compare(new, ref):
    a = pd.read_csv(new, sep="\t", dtype=str).fillna("")
    b = pd.read_csv(ref, sep="\t", dtype=str).fillna("")
    if list(a.columns) != list(b.columns) or a.shape != b.shape:
        return f"열/행 구조가 다름 {a.shape} vs {b.shape}"
    bad = 0
    for c in a.columns:
        x = pd.to_numeric(a[c], errors="coerce")
        y = pd.to_numeric(b[c], errors="coerce")
        num = x.notna() & y.notna()
        if not np.allclose(x[num], y[num], rtol=RTOL, atol=ATOL):
            bad += int((~np.isclose(x[num], y[num], rtol=RTOL, atol=ATOL)).sum())
        bad += int((a.loc[~num, c] != b.loc[~num, c]).sum())
    return None if bad == 0 else f"다른 칸 {bad}개"


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=ROOT / "results")
    ap.add_argument("--update", action="store_true")
    a = ap.parse_args()
    tables = Path(a.results) / "tables"
    ref = Path(a.results) / "reference"
    if not any(tables.glob("*.tsv")):
        sys.exit("results/tables 가 비어 있습니다. 먼저 python run_all.py 를 실행하세요")
    if a.update:
        ref.mkdir(parents=True, exist_ok=True)
        for p in tables.glob("*.tsv"):
            shutil.copy2(p, ref / p.name)
        print(f"기준 결과 저장: {ref}")
        return
    if not any(ref.glob("*.tsv")):
        sys.exit("기준 결과(results/reference)가 없습니다. 실제 데이터로 처음 돌린 뒤 --update 로 저장하세요")
    fails = 0
    for r in sorted(ref.glob("*.tsv")):
        n = tables / r.name
        msg = "결과 파일 없음" if not n.exists() else compare(n, r)
        print(f"  {'일치' if msg is None else '불일치'}  {r.name}" + ("" if msg is None else f"  ({msg})"))
        fails += msg is not None
    if fails:
        sys.exit(f"재현 실패: {fails}개 표가 기준과 다릅니다")
    print("재현 성공: 모든 결과표가 기준과 일치합니다")


if __name__ == "__main__":
    main()
