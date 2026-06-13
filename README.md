# Fast and Scalable Symbolic Maximum Likelihood Estimation for Massive Univariate Histogram Data

[![Status: Under Review](https://img.shields.io/badge/Status-Under_Review-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the data and computational code necessary to fully replicate the analyses, simulations, and figures presented in the working manuscript: **"Fast and Scalable Symbolic Maximum Likelihood Estimation for Massive Univariate Histogram Data"** by Ahmad Hakiim Jamaluddin, Syaiful Anam, Muhammad Jaffri Mohd Nasir, and Andrea Tri Rian Dani.

*Note: A permanent DOI will be minted via Zenodo upon the formal acceptance and publication of the manuscript.*

## Overview

Evaluating classical maximum likelihood estimation (MLE) objectives repeatedly becomes a massive computational bottleneck for ultra-large datasets ($n \ge 10^6$) or within sliding-window architectures on continuous data streams. This project implements a **Symbolic Maximum Likelihood Estimation (SMLE)** framework for univariate continuous data to resolve this limitation. By mapping raw micro-data into histogram-valued symbols bounded by $B$ boundary points ($B-1$ bins), the log-likelihood evaluation cost drops from O(n) to O(B) per optimization step, effectively decoupling estimation cost from sample size.

The codebase evaluates and supports the following components:
* **Asymptotic Infill Benchmarks:** Multi-scenario grids simulating five classic continuous families: Exponential, Normal, Gamma, Weibull, and Lognormal models.
* **Binning Strategy Variations:** Parallel evaluation of Equal-Width (EW) and Equal-Probability (QU) empirical quantile partitioning methods.
* **Streaming Optimization:** An online count update mechanism operating at a constant time cost (O(1) shifts) over sliding windows.
* **Real-Data Infrastructure:** Production-grade validation loops evaluating physiological streams from the MIMIC-III Intensive Care Heart-Rate Feed and the BIG IDEAs Continuous Glucose Monitoring Feed.

## Repository Structure

* `data/`: Data dropzone where raw data assets must be deposited before pipeline execution (contains `mimic_hr.csv` and `bigideas_cgm.csv`).
* `tables/`: Automatically generated directory where the scripts output the publication-ready summary sheets (.csv) including accuracy, runtime scaling matrices, and convergence summaries.
* `figures/`: Automatically generated directory where the scripts output high-DPI vector plots (.pdf) tracking RMSE tracking curves, log-scale scaling curves, and streaming trajectory loops.
* `results/`: Contains localized subdirectory trees (`simulation/` and `realdata/`) storing raw multi-run experimental execution logs.
* `master_replication.R`: The core master pipeline file containing all modular execution steps, optimization hooks, and automated figure rendering.

## System Requirements and Dependencies

All analyses, simulations, and visualizations were programmed and executed in the **R statistical computing environment**. To preserve typeface embedding and exact vector geometry compliance demanded by leading computational statistics journals, visual outputs use a native `cairo_pdf` system driver hook.

To run the master replication script, ensure you have the following packages installed.

**Engine & Aggregation:**
* `data.table` (v1.14+)

**Wrangling Core:**
* `dplyr`
* `tidyr`
* `purrr`
* `tibble`
* `stringr`

**Plotting & Typography:**
* `ggplot2`
* `scales`

## Usage: Replicating the Analysis

To replicate the findings from the paper:

1. Clone this repository to your local machine:
```bash
git clone [https://github.com/HakiimJ/fast-scalable-symbolic-mle.git](https://github.com/HakiimJ/fast-scalable-symbolic-mle.git)
cd fast-scalable-symbolic-mle

```

2. Position data anchors:
Place the filtered empirical asset files (`mimic_hr.csv` and `bigideas_cgm.csv`) inside the `data/` folder.
3. Configure local directory paths:
Open `master_replication.R` in your IDE and point the `mainDir` pointer variable to your exact workspace absolute path:

```R
mainDir <- "/your/absolute/local/path/fast-scalable-symbolic-mle"

```

4. Choose execution configuration:
The pipeline exposes two pre-configured parameter profiles in Section 1. Switch the active configuration to target your execution profile:

```R
# For a quick verification check (R = 10, smaller sample grids)
ACTIVE_CONFIG <- FAST_CONFIG

# For the complete high-fidelity setup used in the text (R = 1000, full matrices)
ACTIVE_CONFIG <- PAPER_CONFIG

```

5. Fire the controller pipeline:
Launch the master execution engine file directly from your terminal console:

```bash
Rscript master_replication.R

```

A reactive status progress monitor tracks real-time generation loop checkpoints directly within the console standard output. Once complete, your updated findings will be formatted into the `tables/` and `figures/` directories.

```
