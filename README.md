# STUDY-21-001175 — Colchicine CTR Statistical Outputs

Public companion site for the Colchicine pilot trial (Protocol 2101175).

## Live pages

| Page | URL |
|------|-----|
| **Statistical results (methods, tables, figures)** | https://alexanderb95-jpg.github.io/colchicine-ctr-stats/ |
| **Statistical code summary** | https://alexanderb95-jpg.github.io/colchicine-ctr-stats/statistical-analysis.html |
| **R script (download, full code)** | https://github.com/alexanderb95-jpg/colchicine-ctr-stats/blob/main/colchicine_ctr_statistical_code.R |
| **Source R Markdown** | https://github.com/alexanderb95-jpg/colchicine-ctr-stats/blob/main/CTR_Submission_Results_Protocol2101175.Rmd |

## Render modes

One R Markdown source (`CTR_Submission_Results_Protocol2101175.Rmd`) supports two outputs via `params$public_site`:

| Mode | Command | Output | Contents |
|------|---------|--------|----------|
| **Public GitHub Pages** | `Rscript scripts/render_public_site.R` | `index.html` | Methods, accrual, characteristics, primary/safety/PK-PD/survival results, figures, statistical analysis, appendix dataset |
| **Full CTR manuscript** | `rmarkdown::render("CTR_Submission_Results_Protocol2101175.Rmd", output_file = "CTR_Submission_Results_Protocol2101175.html")` | `CTR_Submission_Results_Protocol2101175.html` | Public content plus abstract, lessons learned, trial/drug info, discussion, acknowledgments, references |

Narrative sections excluded from the public site live under `sections/` and are included only when `public_site = FALSE`.

## Repository contents

- `index.html` — slim public stats companion (GitHub Pages root)
- `CTR_Submission_Results_Protocol2101175.Rmd` — source R Markdown for all tables/figures
- `sections/` — narrative CTR sections (abstract, discussion, back matter, etc.)
- `scripts/render_public_site.R` — one-command public-site knit
- `scripts/protocol2101175_redcap_helpers.R` — REDCap helper functions
- `colchicine_ctr_statistical_code.R` — concatenated helpers + all R chunks
- `scripts/extract_stats_code.py` — rebuild `colchicine_ctr_statistical_code.R` from the Rmd

## Trial context

- Protocol STUDY-21-001175 (IND 156493), v3.0 (22 May 2023)
- NCT05279690
- Cohort 1 low dose: colchicine 0.6 mg orally twice daily × 14 days
- Primary endpoint: one-sided one-sample *t*-test at α = 0.025 (Bonferroni-adjusted per protocol §12.4)

Place the REDCap export at `data/clinical/protocol_2101175/ProtocolNo2101175PRM_DATA_2026-05-18_2135.csv` before knitting. Patient-level data are **not** committed to this repository.

## Rebuild public site

```bash
# Requires local REDCap CSV (see path above)
Rscript scripts/render_public_site.R
python3 scripts/extract_stats_code.py
git add index.html CTR_Submission_Results_Protocol2101175.Rmd sections/ scripts/render_public_site.R colchicine_ctr_statistical_code.R README.md
git commit -m "feat: slim public stats site with explicit H0 and Bonferroni alpha"
git push
```

## Rebuild full CTR HTML (internal)

```r
rmarkdown::render(
  "CTR_Submission_Results_Protocol2101175.Rmd",
  output_file = "CTR_Submission_Results_Protocol2101175.html",
  params = list(public_site = FALSE)
)
```
