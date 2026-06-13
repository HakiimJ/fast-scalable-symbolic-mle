# Fast and Scalable Symbolic Maximum Likelihood Estimation for Massive Univariate Histogram Data

This repository contains the full data replication pipeline and computational source code necessary to reproduce the analyses, simulations, tables, and figures presented in the manuscript: **"Fast and Scalable Symbolic Maximum Likelihood Estimation for Massive Univariate Histogram Data"** by Ahmad Hakiim Jamaluddin, Syaiful Anam, Muhammad Jaffri Mohd Nasir, and Andrea Tri Rian Dani.

*Note: A permanent DOI will be minted via Zenodo upon the formal acceptance and publication of the manuscript.*

---

## Overview

Maximum likelihood estimation (MLE) is a foundational tool for parametric inference, but evaluating raw likelihoods repeatedly becomes computationally prohibitive for massive samples ($n \ge 10^6$) and sliding-window workflows on streaming data.

This project develops a **Symbolic Maximum Likelihood Estimation (SMLE)** framework for univariate continuous data. By compressing the raw sample into a histogram-valued symbol defined by $B$ boundary points ($B-1$ bins), the evaluation cost of the resulting multinomial log-likelihood is decoupled from $n$ and bounded at $\mathcal{O}(B)$ per optimization step.

### Key Project Components:

* **Asymptotic Approximation Proofs:** Mathematical validation showing that the symbolic log-likelihood converges uniformly to the raw-data log-likelihood up to a parameter-independent bin-width additive term as the partition becomes finer.
* **Multi-Scenario Simulation Study:** Evaluation across five continuous distribution families (Exponential, Normal, Gamma, Weibull, and Lognormal) for sample sizes up to $10^6$ and bin configurations ranging from 10 to 500 bins.
* **Binning Strategy Variations:** Direct comparison between Equal-Width (EW) and Equal-Probability (QU) partitions.
* **Real-Data Streaming Validation:** A sliding-window optimization pipeline with $\mathcal{O}(1)$ streaming count updates applied to the MIMIC-III intensive-care heart-rate stream (containing 1,912,844 observations).

---

## Repository Structure

```text
fast-scalable-symbolic-mle/
├── master_replication.R       # The master script containing the full pipeline
├── data/                      # Local data directory for real-world validation
│   └── mimic_hr.csv           # Prepared MIMIC-III heart rate stream data
├── figures/                   # Automatically generated manuscript-ready PDF plots
│   ├── Fig_Sim_RMSE_vs_B.pdf
│   ├── Fig_Sim_Runtime_vs_n.pdf
│   └── Fig_MIMIC_Window_Params.pdf
├── tables/                    # Automatically generated manuscript-ready CSV tables
│   ├── Tab_Sim_Accuracy_1.csv
│   ├── Tab_Sim_Runtime_1.csv
│   └── Tab_RealData_Runtime.csv
└── results/                   # Detailed replication artifacts and checkpoint data
    ├── simulation/
    └── realdata/

```

---

## System Requirements and Dependencies

All analyses, benchmarks, and data visualisations were programmed and executed within the **R statistical computing environment** (v4.0 or later recommended).

The pipeline includes an automated dependency manager (`load_or_install`) that checks your local environment and fetches missing packages from CRAN. The script relies on the following packages:

* **Data Wrangling & I/O:** `data.table`, `dplyr`, `tidyr`, `purrr`, `tibble`, `stringr`
* **Visualization & Formatting:** `ggplot2`, `scales`
* **Graphics Device Extension:** Output figures are rendered via `cairo_pdf` to ensure vector-graphic manuscript compliance.

---

## Usage: Replicating the Analysis

To replicate all simulation tables, NLL gap assessments, convergence frequencies, and real-data window trajectories from the paper:

1. **Clone the repository:**
```bash

```



git clone https://github.com/HakiimJ/fast-scalable-symbolic-mle.git
cd fast-scalable-symbolic-mle

```

2. **Configure Local Environment Paths:**
   Open `master_replication.R` in your preferred editor/IDE and change the `mainDir` path at the top of the script to reflect your local directory structure:
   ```R
mainDir <- "/your/local/path/to/fast-scalable-symbolic-mle"

```

3. **Configure Runtime Modes (Optional):**
By default, the script runs the complete `PAPER_CONFIG` used in the final manuscript ($R = 1,000$ replications, large sample matrices, full window evaluation). If you want to check your system capabilities or verify execution quickly, swap the active config flag to `FAST_CONFIG`:
```R

```



ACTIVE_CONFIG <- FAST_CONFIG  # Short test run

```

4. **Execute the full pipeline:**
   Run the master controller from your shell terminal:
   ```bash
Rscript master_replication.R

```

Or execute `run_full_pipeline(ACTIVE_CONFIG)` inside an interactive R session/RStudio.

The console will report progress percentages with a live tracking bar. Upon completion, all manuscript-ready tables and figures will be systematically exported into their designated folders.
