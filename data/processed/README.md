# 분석용 데이터 (GitHub에 올림)

`scripts/prepare_data.py` 가 만듭니다. 한 번 만들어 커밋해 두면, 이후에는 어느 컴퓨터에서든
인터넷 없이 `python run_all.py` 만으로 분석이 재현됩니다.

| 파일 | 내용 |
|---|---|
| `tpm.tsv` | 유전자 기호 × 15개 샘플 TPM |
| `counts.tsv` | 유전자 기호 × 15개 샘플 expected count |
| `samples.tsv` | sample, individual(A–E), stage(W/I/B) |
| `prepare_report.txt` | 처리 기록 (파일 수, 주석 방법, 유전자 수) |

원자료: Ferreira et al. 2017, Dryad doi:10.5061/dryad.7v1m1
