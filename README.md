# osip-mini-project

**How does noise in arterial blood data and PET time-activity curves propagate into kinetic-model estimates of the total volume of distribution (V<sub>T</sub>)?**

A small, fully reproducible simulation study. All data are simulated; no real subject or study data are used.

Mini-project for the KI PhD course *Open Science in Practice: Collaborative Research with Git* (K8F6106, HT 2026).

> **Status:** planning. The sections below describe the intended pipeline; they will be updated as the code is written.

## Question

For a reversibly binding radioligand, how large are the bias and variability of V<sub>T</sub> when:

1. only the **blood/plasma input function** is noisy,
2. only the **tissue TAC** is noisy,
3. both are noisy?

And how does this differ between a **two-tissue compartment model (2TCM)** and **Logan graphical analysis**?

## Pipeline

| Step | What it does | Output |
|---|---|---|
| 1. Input function | Generate a noise-free arterial plasma curve from a parametric model (e.g. Feng), sampled at a realistic blood-sampling schedule | `data/processed/` |
| 2. Tissue TACs | Simulate noise-free TACs with a 2TCM for a set of known rate constants (K<sub>1</sub>–k<sub>4</sub>), integrated over a realistic PET frame schedule | `data/processed/` |
| 3. Noise | Add noise to the blood samples (measurement error, sampling-time error) and/or to the TACs (frame-duration- and count-dependent), at several noise levels, with many replicates per level | `data/processed/` |
| 4. Fitting | Estimate V<sub>T</sub> from each noisy replicate with 2TCM and Logan | `results/tables/` |
| 5. Summary | Bias and coefficient of variation of V<sub>T</sub> per noise scenario and method | `results/tables/`, `results/figures/` |

Planned main result: a figure of V<sub>T</sub> bias and CV against noise level for each scenario and method, plus a summary table.

## Design choices (to be fixed in code)

- Simulation parameters (rate constants, frame schedule, blood sampling times, noise levels, number of replicates) are kept in **one parameter file**, not scattered across scripts.
- All random steps use a **fixed seed**, so the full pipeline gives identical results on every run.
- Each pipeline step is a **separate script**, run in numbered order; intermediate data are written to `data/processed/` and can be regenerated at any time.
- No step is done by hand; everything in `results/` is produced by the scripts.

## Structure

```
data/raw/          (unused — all data are simulated)
data/processed/    simulated input functions, TACs and noisy replicates
scripts/           pipeline scripts, run in numbered order
results/figures/   figures produced by the pipeline
results/tables/    tables produced by the pipeline
```

## How to reproduce

_To be written once the language and environment setup are decided._
