# Fast and Scalable Symbolic Maximum Likelihood Estimation for Massive Univariate Histogram Data

This repository hosts the replication code and data pipeline for the manuscript: **"Fast and Scalable Symbolic Maximum Likelihood Estimation for Massive Univariate Histogram Data"** authored by Ahmad Hakiim Jamaluddin, Syaiful Anam, Muhammad Jaffri Mohd Nasir, and Andrea Tri Rian Dani.

---

## Overview

Evaluating classical maximum likelihood estimation (MLE) objectives repeatedly becomes a massive bottleneck for ultra-large datasets ($n \ge 10^6$) or within sliding-window architectures on continuous data streams.

This project implements a **Symbolic Maximum Likelihood Estimation (SMLE)** framework for univariate continuous data. By mapping raw micro-data into histogram-valued symbols bounded by $B$ boundary points ($B-1$ bins), the log-likelihood evaluation cost drops from $\mathcal{O}(n)$ to $\mathcal{O}(B)$ per optimization step.

### What This Pipeline Does:

* **Asymptotic Infill Benchmarks:** Simulates five classic continuous families (Exponential, Normal, Gamma, Weibull, Lognormal) up to $n=10^6$ under both Equal-Width (EW) and Equal-Probability (QU) binning schemes.
* **Real-Data Moving-Window Applications:** Implements an optimized online count update mechanism ($\mathcal{O}(1)$ count shifts) evaluated across two physiological monitoring feeds:
1. **MIMIC-III Intensive Care Heart-Rate Feed** (*The primary empirical case presented in the manuscript text*).
2. **BIG IDEAs Continuous Glucose Monitoring Feed** (*Maintained within the source code pipeline execution structure*).



---

## Repository Architecture

The execution of the replication script automatically generates a clean directory tree to separate structural components, intermediate outputs, and publication-ready outputs:

```text
fast-scalable-symbolic-mle/
├── master_replication.R       # Single cohesive script running the full workflow
├── data/                      # Data dropzone (Populate before execution)
│   ├── mimic_hr.csv           # Filtered heart-rate numerical stream data
│   └── bigideas_cgm.csv       # Filtered glucose numerical stream data
├── tables/                    # Publication-ready summary sheets (.csv)
│   ├── Tab_Sim_Accuracy_1.csv
│   ├── Tab_Sim_Accuracy_2.csv
│   ├── Tab_Sim_Runtime_1.csv
│   ├── Tab_Sim_Runtime_2.csv
│   ├── Tab_Sim_Convergence_and_NLLGap.csv
│   ├── Tab_RealData_Runtime.csv
│   └── Tab_Manuscript_Summary_Stats.csv
├── figures/                   # High-DPI standalone vector plots (.pdf)
│   ├── Fig_Sim_RMSE_vs_B.pdf
│   ├── Fig_Sim_Runtime_vs_n.pdf
│   ├── Fig_MIMIC_Window_Params.pdf
│   ├── Fig_BIGIDEAs_Window_Params.pdf
│   └── Fig_RealData_Sensitivity.pdf
└── results/                   # Raw multi-run experimental logs
    ├── simulation/
    └── realdata/

```

---

## Prerequisites & Setup

The entire computational suite is engineered in plain **R**. It includes a lightweight bootstrapping function (`load_or_install`) that dynamically flags, downloads, and maps dependencies directly from the cloud CRAN mirror without requiring manual package management.

### Required Packages:

* **Engine & Aggregation:** `data.table`
* **Wrangling Core:** `dplyr`, `tidyr`, `purrr`, `tibble`, `stringr`
* **Plotting and Typography:** `ggplot2`, `scales`

*Note: Visual outputs use a `cairo_pdf` device hook to guarantee embedded typefaces and strict vector geometry compliance required by leading computational statistics journals.*

---

## Replicating the Results

### 1. Position Data Anchors

Drop your filtered empirical asset sheets (`mimic_hr.csv` and `bigideas_cgm.csv`) inside the newly established `data/` folder.

### 2. Configure Local Directory Paths

Open `master_replication.R` and point the `mainDir` pointer variable to your exact workspace absolute path:

```R
mainDir <- "/your/absolute/local/path/fast-scalable-symbolic-mle"

```

### 3. Choose Execution Configuration

The script exposes two preset configuration matrices under Section 1:

* `FAST_CONFIG`: Tailored for local validation loops ($R = 10$, truncated samples, micro-window depth).
* `PAPER_CONFIG`: The complete, high-fidelity setup utilized in the finalized text ($R = 1000$ iterations, full sample matrices, absolute streaming depth).

Switch the active assignment flag to target your resource profile:

```R
ACTIVE_CONFIG <- PAPER_CONFIG  # Run final paper specs

```

### 4. Fire the Controller Pipeline

Launch the master file directly from your terminal console:

```bash
Rscript master_replication.R

```

A reactive status monitor tracks real-time generation metrics directly within the console standard output. Once finalized, look in the `tables/` and `figures/` directories to check the compiled outputs.
