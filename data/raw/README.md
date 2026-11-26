# 원자료 (GitHub에는 올리지 않음)

Dryad: https://doi.org/10.5061/dryad.7v1m1 (Ferreira et al. 2017, Molecular Ecology)

이 폴더에 아래 파일을 브라우저로 내려받아 두세요.

| 파일 | 크기 | 필수 |
|---|---|---|
| `Lep_ame_RSEM.out.genes.results.tar` | 88 MB | 필수 |
| `README_for_Lep_ame_RSEM.out.genes.results.rtf` | 1 KB | 권장 (열 설명 확인) |
| `Lep_ame_filtered_skin_transcriptome.fasta` | 190 MB | gene_id 가 Trinity ID 일 때 |
| `biomart_rabbit.tsv`, `biomart_mouse.tsv` | - | Ensembl ID → 유전자 기호 변환 |

BioMart 표 만드는 법: https://www.ensembl.org/biomart → Dataset 에서 *Rabbit genes* (다음에 *Mouse genes* 로 한 번 더)
→ Attributes 에서 `Gene stable ID`, `Transcript stable ID`, `Protein stable ID`, `Gene name` 선택 → Results → TSV 로 내보내기.
