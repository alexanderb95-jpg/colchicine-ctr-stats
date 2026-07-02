# STUDY-21-001175 — Colchicine CTR Statistical Outputs

Public companion site for the Colchicine pilot trial CTR submission.

## Live pages

| Page | URL |
|------|-----|
| **Full trial outputs (tables + figures)** | https://alexanderb95-jpg.github.io/colchicine-ctr-stats/ |
| **Statistical code summary** | https://alexanderb95-jpg.github.io/colchicine-ctr-stats/statistical-analysis.html |
| **R script (download, full code)** | https://github.com/alexanderb95-jpg/colchicine-ctr-stats/blob/main/colchicine_ctr_statistical_code.R |
| **Source R Markdown** | https://github.com/alexanderb95-jpg/colchicine-ctr-stats/blob/main/CTR_Submission_Results_Protocol2101175.Rmd |

## Repository contents

- `index.html` — PI-verified outcomes report + complete statistical code appendix at bottom
- `CTR_Submission_Results_Protocol2101175.Rmd` — source R Markdown for all tables/figures
- `protocol2101175_redcap_helpers.R` — REDCap helper functions (screen failure, survival, CTCAE labels)
- `colchicine_ctr_statistical_code.R` — concatenated helpers + all R chunks (~1,660 lines)
- `scripts/extract_stats_code.py` — rebuild `colchicine_ctr_statistical_code.R` from the Rmd

## Trial context

- Protocol STUDY-21-001175 (IND 156493), v3.0 (22 May 2023)
- NCT05279690
- Cohort 1 low dose: colchicine 0.6 mg orally twice daily × 14 days
- Data cut: 23 May 2026

No patient-level export data are included in this repository.

## Rebuild statistical code file

```bash
python3 scripts/extract_stats_code.py
git add colchicine_ctr_statistical_code.R index.html statistical-analysis.html
git commit -m "chore: refresh extracted statistical code"
git push
```
