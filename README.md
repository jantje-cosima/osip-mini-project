# Receptor occupancy and plasma EC50 – simulated PET blocking study

Mini-project for the KI PhD course *Open Science in Practice: Collaborative
Research with Git* (K8F6106, HT 2026).

A reproducible R pipeline that estimates receptor occupancy with the **Lassen
plot** and the plasma **EC50** with an Emax model, with quality control (QC)
after every step. All data are **simulated**; no real subject data are used.

![EC50 curve](results/figures/ec50_curve.png)

## Data

| File | Content |
|---|---|
| `data/raw/vt_estimates.csv` | VT (2TCM) per subject × session × region, with SE and fit flags |
| `data/raw/drug_plasma.csv` | plasma drug concentration at the start and end of each block scan |
| `data/raw/tacs.csv` | regional time-activity curves (optional QC) |
| `data/raw/plasma_input.csv` | metabolite-corrected arterial input (optional QC) |

10 subjects × 3 sessions (baseline, block 1 at ~3 h, block 2 at ~24 h after a
single oral dose of 1–100 mg) × 10 brain regions.

## Pipeline

```
scripts/
├── 00_run_pipeline.R      run_pipeline() – runs every stage and renders the reports
├── 01_load_data.R         load_vt(), load_drug()      – read + structural checks
├── 02_qc_input.R          qc_vt(), qc_drug()          – QC of input data
├── 03_drug_exposure.R     drug_exposure()             – plasma concentration per scan
├── 04_fit_lassen.R        fit_lassen()                – occupancy + VND per scan
├── 05_qc_lassen.R         qc_lassen()                 – QC of Lassen fits
├── 06_fit_ec50.R          fit_ec50()                  – Emax model
├── 07_qc_ec50.R           qc_ec50()                   – QC of EC50 fit
├── 08_make_outputs.R      make_outputs()              – figures and tables
├── qc_utils.R             shared QC helpers (flag format)
└── optional_refit_vt.R    load + QC of TACs and arterial input (optional)

reports/
├── report.qmd             final report
└── qc_report.qmd          QC report (all flags, all Lassen fits)

results/
├── figures/               final figures
└── tables/                final tables (occupancy per scan, EC50, QC flags)
```

A flowchart of the pipeline is in `docs/pipeline.drawio`.

QC never deletes data: each check flags rows as *exclude* (not used in later
steps) or *warning* (kept, reported). All flags are collected in the QC report.

## How to reproduce

Requirements: **R ≥ 4.6** and **Quarto** (bundled with Positron, or from
<https://quarto.org>).

```bash
git clone git@github.com:marbol1/osip-mini-project.git
cd osip-mini-project
```

In R, from the project folder, install the exact package versions:

```r
renv::restore()
```

Then run the whole pipeline (from the terminal):

```bash
Rscript scripts/00_run_pipeline.R
```

or interactively in R:

```r
source("scripts/00_run_pipeline.R")
res <- run_pipeline()                 # main pipeline + both reports
res <- run_pipeline(tac_qc = TRUE)    # also QC of TACs and arterial input
```

This writes `results/figures/`, `results/tables/`, `reports/report.html` and
`reports/qc_report.html`.

## Results

Plasma EC50 = **23.0 ng/mL** (95 % CI 19.3–27.3), from 20 block scans
(see `results/tables/ec50_result.csv` and the final report).