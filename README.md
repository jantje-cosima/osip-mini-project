# osip-mini-project

Mini-project for the KI PhD course *Open Science in Practice: Collaborative Research with Git* (K8F6106, HT 2026).

## Structure

```
data/raw/          input data (never edited by hand)
data/processed/    intermediate data produced by the pipeline
scripts/           analysis scripts, run in numbered order

scripts/
├── 00_run_pipeline.R      run_pipeline()   – master: sources all scripts, runs every stage, renders the report
├── 01_load_data.R         load_vt(), load_drug()                – read + check columns/types/units
├── 02_qc_input.R          qc_vt(), qc_drug()                    – QC of input data (incl. VT fit quality)
├── 03_drug_exposure.R     drug_exposure()                       – mean plasma drug conc. per scan
├── 04_fit_lassen.R        fit_lassen()                          – occupancy + VND per subject × block
├── 05_qc_lassen.R         qc_lassen()                           – QC of Lassen fits
├── 06_fit_ec50.R          fit_ec50()                            – Emax model: occupancy vs plasma conc.
├── 07_qc_ec50.R           qc_ec50()                             – QC of EC50 fit
├── 08_make_outputs.R      figures + tables to results/
├── qc_utils.R             helpers shared by all qc_*() – collect flags into the QC report
└── optional_refit_vt.R    load_tacs(), load_plasma(), qc_tacs(), qc_plasma(), refit_vt()
                           – optional 2TCM refit (run_pipeline(refit = TRUE))

report.qmd                 final report (project root, rendered by run_pipeline())

results/figures/   figures produced by the pipeline
results/tables/    tables/results produced by the pipeline
```

## How to reproduce

_To be written._
