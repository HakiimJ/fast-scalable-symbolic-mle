options(stringsAsFactors = FALSE)

# ============================================================================
# FULL REPLICATION SCRIPT
# "Fast and Scalable Symbolic Maximum Likelihood Estimation for Massive Univariate Histogram Data"
# ============================================================================
# PURPOSE:
#   - Reproduce simulation tables and figures for the manuscript
#   - Run sliding-window real-data benchmarks on strict local CSV datasets
#   - Export manuscript-ready CSV tables and PDF figures to disk
#
# NOTE:
#   Update `mainDir` below to your local manuscript directory before running.
#   Make sure to download the real dataset and modify the directory code accordingly
# ============================================================================

# ==============================
# 0. USER PATHS AND SETUP
# ==============================
mainDir <- getwd()

for (sub in c("", "figures", "tables", "data", "results", "results/simulation", "results/realdata")) {
  dir.create(file.path(mainDir, sub), showWarnings = FALSE, recursive = TRUE)
}

load_or_install <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "scales",
  "purrr", "tibble", "stringr"
)
load_or_install(required_pkgs)

set.seed(20260101)

cat("\n====================================================================\n")
cat(" SYMBOLIC MAXIMUM LIKELIHOOD ESTIMATION: FULL REPLICATION SCRIPT\n")
cat(sprintf(" Root Path: %s\n", mainDir))
cat("====================================================================\n\n")

# ==============================
# 1. CONFIGURATION
# ==============================
FAST_CONFIG <- list(
  n_grid           = c(1e4, 1e5),
  B_grid           = c(10, 20, 50, 100),
  reps             = 10,
  families         = c("exp", "normal", "gamma", "weibull", "lognormal"),
  stream_win_hr    = 120L,
  stream_step_hr   = 10L,
  stream_win_cgm   = 48L,
  stream_step_cgm  = 5L,
  sens_B_grid      = c(10, 20, 50),
  sens_W_grid_hr   = c(60, 120),
  sens_W_grid_cgm  = c(24, 48),
  stream_max_win   = 50L,
  realdata_B       = 20L,
  run_full_windows = FALSE
)

PAPER_CONFIG <- list(
  n_grid           = c(1e4, 1e5, 1e6),
  B_grid           = c(10, 20, 50, 100, 200, 500),
  reps             = 1000,
  families         = c("exp", "normal", "gamma", "weibull", "lognormal"),
  stream_win_hr    = 10000L, # Increased to reveal symbolic speed advantage
  stream_step_hr   = 500L,
  stream_win_cgm   = 1000L,  # Increased for BIG IDEAs (Total N is ~2500)
  stream_step_cgm  = 50L,
  sens_B_grid      = c(5, 10, 15, 20),
  sens_W_grid_hr   = c(5000, 10000, 20000),
  sens_W_grid_cgm  = c(500, 1000, 1500),
  stream_max_win   = Inf,
  realdata_B       = 20L,
  run_full_windows = TRUE
)

# Changed to PAPER_CONFIG for final run
ACTIVE_CONFIG <- PAPER_CONFIG

EPS_PROB <- 1e-12
OPT_MAXIT <- 1000

# ==============================
# 2. GLOBAL UTILITIES
# ==============================
fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "f"))
}

safe_divide <- function(num, den) {
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

save_plot <- function(plot_obj, filename, width = 10, height = 6) {
  ggsave(
    filename = file.path(mainDir, "figures", filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300,
    device = cairo_pdf
  )
}

safe_write_csv <- function(x, filename, folder = "tables") {
  write.csv(x, file.path(mainDir, folder, filename), row.names = FALSE)
}

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold", size = rel(1.1), margin = margin(b = 6), hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = rel(0.9))
    )
}

model_pretty <- c(
  exp = "Exponential",
  normal = "Normal",
  gamma = "Gamma",
  weibull = "Weibull",
  lognormal = "Lognormal"
)

strategy_pretty <- c(
  equal_width = "Equal Width",
  quantile = "Equal Probability"
)

# ==============================
# 3. DATA GENERATION FOR SIMULATION
# ==============================
get_truth <- function(model) {
  switch(
    model,
    exp       = c(rate = 1.2),
    normal    = c(mu = 0, sigma = 1),
    gamma     = c(shape = 2, rate = 1),
    weibull   = c(shape = 1.5, scale = 1),
    lognormal = c(meanlog = 0, sdlog = 0.75),
    stop("Unknown model: ", model)
  )
}

generate_data <- function(n, model) {
  switch(
    model,
    exp       = rexp(n, rate = 1.2),
    normal    = rnorm(n, mean = 0, sd = 1),
    gamma     = rgamma(n, shape = 2, rate = 1),
    weibull   = rweibull(n, shape = 1.5, scale = 1),
    lognormal = rlnorm(n, meanlog = 0, sdlog = 0.75),
    stop("Unknown model: ", model)
  )
}

# ==============================
# 4. HISTOGRAM CONSTRUCTION
# ==============================
make_histogram <- function(y, B, strategy = c("quantile", "equal_width")) {
  strategy <- match.arg(strategy)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  stopifnot(n > 0, B >= 2)
  
  if (strategy == "equal_width") {
    rng <- range(y)
    if (rng[1] == rng[2]) {
      rng[1] <- rng[1] - 0.5
      rng[2] <- rng[2] + 0.5
    }
    boundaries <- seq(rng[1], rng[2], length.out = B + 1)
  } else {
    probs <- seq(0, 1, length.out = B + 1)
    boundaries <- unname(quantile(y, probs = probs, type = 8, names = FALSE, na.rm = TRUE))
    if (anyDuplicated(boundaries)) {
      for (i in 2:length(boundaries)) {
        if (boundaries[i] <= boundaries[i - 1]) boundaries[i] <- boundaries[i - 1] + 1e-8
      }
    }
  }
  
  boundaries[1] <- boundaries[1] - 1e-8
  boundaries[length(boundaries)] <- boundaries[length(boundaries)] + 1e-8
  
  breaks_factor <- cut(y, breaks = boundaries, include.lowest = TRUE, right = FALSE)
  counts <- as.integer(tabulate(as.integer(breaks_factor), nbins = B))
  
  list(
    boundaries = boundaries,
    counts = counts,
    n = n,
    B = B,
    strategy = strategy,
    widths = diff(boundaries)
  )
}

find_bin_index <- function(value, boundaries) {
  idx <- findInterval(value, boundaries, rightmost.closed = FALSE, all.inside = TRUE)
  idx[idx < 1] <- 1L
  idx[idx >= length(boundaries)] <- length(boundaries) - 1L
  as.integer(idx)
}

# ==============================
# 5. LIKELIHOODS AND FITTING
# ==============================
transform_params <- function(z, model) {
  switch(
    model,
    normal    = c(mu = z[1], sigma = exp(z[2])),
    exp       = c(rate = exp(z[1])),
    gamma     = c(shape = exp(z[1]), rate = exp(z[2])),
    weibull   = c(shape = exp(z[1]), scale = exp(z[2])),
    lognormal = c(meanlog = z[1], sdlog = exp(z[2])),
    stop("Unknown model: ", model)
  )
}

get_init <- function(y, model) {
  y <- as.numeric(y)
  y_pos <- y[y > 0 & is.finite(y)]
  switch(
    model,
    normal = {
      s <- sd(y)
      c(mean(y), log(ifelse(is.finite(s) && s > 0, s, 1)))
    },
    exp = {
      m <- mean(y[y > 0])
      c(log(1 / max(m, 1e-8)))
    },
    gamma = {
      m <- mean(y_pos)
      v <- var(y_pos)
      v <- ifelse(is.finite(v) && v > 0, v, 1)
      c(log(max(m^2 / v, 1e-8)), log(max(m / v, 1e-8)))
    },
    weibull = c(log(1.1), log(max(mean(y_pos), 1e-8))),
    lognormal = {
      ly <- log(y_pos)
      s <- sd(ly)
      c(mean(ly), log(ifelse(is.finite(s) && s > 0, s, 0.5)))
    },
    stop("Unknown model: ", model)
  )
}

raw_nll <- function(y, model, z) {
  pars <- transform_params(z, model)
  switch(
    model,
    normal    = -sum(dnorm(y, mean = pars["mu"], sd = pars["sigma"], log = TRUE)),
    exp       = -sum(dexp(y, rate = pars["rate"], log = TRUE)),
    gamma     = -sum(dgamma(y, shape = pars["shape"], rate = pars["rate"], log = TRUE)),
    weibull   = -sum(dweibull(y, shape = pars["shape"], scale = pars["scale"], log = TRUE)),
    lognormal = -sum(dlnorm(y, meanlog = pars["meanlog"], sdlog = pars["sdlog"], log = TRUE)),
    stop("Unknown model: ", model)
  )
}

symbolic_probs <- function(boundaries, model, pars) {
  left <- boundaries[-length(boundaries)]
  right <- boundaries[-1]
  probs <- switch(
    model,
    normal    = pnorm(right, mean = pars["mu"], sd = pars["sigma"]) - pnorm(left, mean = pars["mu"], sd = pars["sigma"]),
    exp       = pexp(right, rate = pars["rate"]) - pexp(left, rate = pars["rate"]),
    gamma     = pgamma(right, shape = pars["shape"], rate = pars["rate"]) - pgamma(left, shape = pars["shape"], rate = pars["rate"]),
    weibull   = pweibull(right, shape = pars["shape"], scale = pars["scale"]) - pweibull(left, shape = pars["shape"], scale = pars["scale"]),
    lognormal = plnorm(right, meanlog = pars["meanlog"], sdlog = pars["sdlog"]) - plnorm(left, meanlog = pars["meanlog"], sdlog = pars["sdlog"]),
    stop("Unknown model: ", model)
  )
  pmax(probs, EPS_PROB)
}

symbolic_nll <- function(hist_obj, model, z) {
  pars <- transform_params(z, model)
  probs <- symbolic_probs(hist_obj$boundaries, model, pars)
  -sum(hist_obj$counts * log(probs))
}

fit_mle <- function(y, model, B = NULL, strategy = NULL, init = NULL) {
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (is.null(init)) init <- get_init(y, model)
  
  if (is.null(B)) {
    obj <- function(z) raw_nll(y, model, z)
    t0 <- proc.time()["elapsed"]
    opt <- tryCatch(
      suppressWarnings(optim(init, obj, method = "BFGS", control = list(maxit = OPT_MAXIT))),
      error = function(e) list(par = init, value = NA_real_, convergence = 99, message = conditionMessage(e))
    )
    elapsed <- proc.time()["elapsed"] - t0
    est <- transform_params(opt$par, model)
    return(list(
      est = est,
      nll = opt$value,
      converged = isTRUE(opt$convergence == 0),
      time = as.numeric(elapsed),
      time_hist = 0,
      time_opt = as.numeric(elapsed),
      hist_obj = NULL
    ))
  }
  
  t_hist0 <- proc.time()["elapsed"]
  hist_obj <- make_histogram(y, B = B, strategy = strategy)
  t_hist <- proc.time()["elapsed"] - t_hist0
  
  obj <- function(z) symbolic_nll(hist_obj, model, z)
  t_opt0 <- proc.time()["elapsed"]
  opt <- tryCatch(
    suppressWarnings(optim(init, obj, method = "BFGS", control = list(maxit = OPT_MAXIT))),
    error = function(e) list(par = init, value = NA_real_, convergence = 99, message = conditionMessage(e))
  )
  t_opt <- proc.time()["elapsed"] - t_opt0
  
  est <- transform_params(opt$par, model)
  list(
    est = est,
    nll = opt$value,
    converged = isTRUE(opt$convergence == 0),
    time = as.numeric(t_hist + t_opt),
    time_hist = as.numeric(t_hist),
    time_opt = as.numeric(t_opt),
    hist_obj = hist_obj
  )
}

param_rmse <- function(est, truth) {
  est <- as.numeric(est)
  truth <- as.numeric(truth)
  sqrt(mean((est - truth)^2, na.rm = TRUE))
}

param_absdiff <- function(est1, est2) {
  est1 <- as.numeric(est1)
  est2 <- as.numeric(est2)
  mean(abs(est1 - est2), na.rm = TRUE)
}

# ==============================
# 6. SIMULATION STUDY
# ==============================
run_one_sim <- function(model, n, B, rep_id) {
  set.seed(202600000 + rep_id + as.integer(n) + match(model, names(model_pretty)) * 1000 + B)
  y <- generate_data(n, model)
  truth <- get_truth(model)
  
  res_raw <- fit_mle(y, model)
  res_ew  <- fit_mle(y, model, B = B, strategy = "equal_width", init = get_init(y, model))
  res_qu  <- fit_mle(y, model, B = B, strategy = "quantile", init = get_init(y, model))
  
  raw_nll_at_raw <- raw_nll(y, model, z = unname(c(
    if (model == "normal") c(res_raw$est["mu"], log(res_raw$est["sigma"])) else NULL,
    if (model == "exp") log(res_raw$est["rate"]) else NULL,
    if (model == "gamma") c(log(res_raw$est["shape"]), log(res_raw$est["rate"])) else NULL,
    if (model == "weibull") c(log(res_raw$est["shape"]), log(res_raw$est["scale"])) else NULL,
    if (model == "lognormal") c(res_raw$est["meanlog"], log(res_raw$est["sdlog"])) else NULL
  )))
  
  raw_nll_at_ew <- raw_nll(y, model, z = unname(c(
    if (model == "normal") c(res_ew$est["mu"], log(res_ew$est["sigma"])) else NULL,
    if (model == "exp") log(res_ew$est["rate"]) else NULL,
    if (model == "gamma") c(log(res_ew$est["shape"]), log(res_ew$est["rate"])) else NULL,
    if (model == "weibull") c(log(res_ew$est["shape"]), log(res_ew$est["scale"])) else NULL,
    if (model == "lognormal") c(res_ew$est["meanlog"], log(res_ew$est["sdlog"])) else NULL
  )))
  
  raw_nll_at_qu <- raw_nll(y, model, z = unname(c(
    if (model == "normal") c(res_qu$est["mu"], log(res_qu$est["sigma"])) else NULL,
    if (model == "exp") log(res_qu$est["rate"]) else NULL,
    if (model == "gamma") c(log(res_qu$est["shape"]), log(res_qu$est["rate"])) else NULL,
    if (model == "weibull") c(log(res_qu$est["shape"]), log(res_qu$est["scale"])) else NULL,
    if (model == "lognormal") c(res_qu$est["meanlog"], log(res_qu$est["sdlog"])) else NULL
  )))
  
  tibble(
    model = model,
    n = n,
    B = B,
    rep = rep_id,
    rmse_raw = param_rmse(res_raw$est, truth),
    rmse_ew = param_rmse(res_ew$est, truth),
    rmse_qu = param_rmse(res_qu$est, truth),
    rel_rmse_ew = safe_divide(param_rmse(res_ew$est, truth), param_rmse(res_raw$est, truth)),
    rel_rmse_qu = safe_divide(param_rmse(res_qu$est, truth), param_rmse(res_raw$est, truth)),
    nll_gap_ew = raw_nll_at_ew - raw_nll_at_raw,
    nll_gap_qu = raw_nll_at_qu - raw_nll_at_raw,
    conv_raw = res_raw$converged,
    conv_ew = res_ew$converged,
    conv_qu = res_qu$converged,
    time_raw = res_raw$time,
    time_ew = res_ew$time,
    time_qu = res_qu$time,
    speedup_ew = safe_divide(res_raw$time, res_ew$time),
    speedup_qu = safe_divide(res_raw$time, res_qu$time)
  )
}

run_simulation_study <- function(cfg) {
  cat("Running simulation study...\n")
  jobs <- expand.grid(
    model = cfg$families,
    n = cfg$n_grid,
    B = cfg$B_grid,
    rep = seq_len(cfg$reps),
    stringsAsFactors = FALSE
  )
  
  out <- vector("list", nrow(jobs))
  total_jobs <- nrow(jobs)
  bar_width <- 50
  
  for (i in seq_len(total_jobs)) {
    job <- jobs[i, ]
    out[[i]] <- run_one_sim(job$model, job$n, job$B, job$rep)
    
    # Progress bar monitor
    pct <- i / total_jobs
    filled <- round(pct * bar_width)
    empty <- bar_width - filled
    bar <- paste0(rep("=", filled), collapse = "")
    if (empty > 0) {
      bar <- paste0(bar, ">")
      empty <- empty - 1
    }
    space <- paste0(rep("-", max(0, empty)), collapse = "")
    
    cat(sprintf("\r  Progress: [%s%s] %d / %d (%5.1f%%)", bar, space, i, total_jobs, 100 * pct))
    flush.console()
  }
  cat("\n") # New line after 100% completion
  
  sim_df <- bind_rows(out)
  safe_write_csv(sim_df, "simulation_replication_all_reps.csv", folder = "results/simulation")
  sim_df
}

generate_sim_tables <- function(sim_df) {
  cat("Generating simulation tables...\n")
  
  tab_acc <- sim_df %>%
    group_by(model, n, B) %>%
    summarise(
      `RMSE (Raw)` = mean(rmse_raw, na.rm = TRUE),
      `RMSE (Sym-EW)` = mean(rmse_ew, na.rm = TRUE),
      `RMSE (Sym-QU)` = mean(rmse_qu, na.rm = TRUE),
      `Rel. RMSE (EW)` = mean(rel_rmse_ew, na.rm = TRUE),
      `Rel. RMSE (QU)` = mean(rel_rmse_qu, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Model = model) %>%
    select(Model, n, B, everything(), -model)
  
  tab_acc_1 <- tab_acc %>% filter(Model %in% c("exp", "lognormal", "normal"))
  tab_acc_2 <- tab_acc %>% filter(Model %in% c("gamma", "weibull"))
  
  safe_write_csv(tab_acc_1, "Tab_Sim_Accuracy_1.csv")
  safe_write_csv(tab_acc_2, "Tab_Sim_Accuracy_2.csv")
  
  tab_rt <- sim_df %>%
    group_by(model, n, B) %>%
    summarise(
      `Time (Raw)` = mean(time_raw, na.rm = TRUE),
      `Time (Sym-EW)` = mean(time_ew, na.rm = TRUE),
      `Time (Sym-QU)` = mean(time_qu, na.rm = TRUE),
      `Speed-up (EW)` = mean(speedup_ew, na.rm = TRUE),
      `Speed-up (QU)` = mean(speedup_qu, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Model = model) %>%
    select(Model, n, B, everything(), -model)
  
  tab_rt_1 <- tab_rt %>% filter(Model %in% c("exp", "lognormal", "normal"))
  tab_rt_2 <- tab_rt %>% filter(Model %in% c("gamma", "weibull"))
  
  safe_write_csv(tab_rt_1, "Tab_Sim_Runtime_1.csv")
  safe_write_csv(tab_rt_2, "Tab_Sim_Runtime_2.csv")
  
  tab_conv <- sim_df %>%
    group_by(model, n, B) %>%
    summarise(
      conv_raw_pct = 100 * mean(conv_raw, na.rm = TRUE),
      conv_ew_pct = 100 * mean(conv_ew, na.rm = TRUE),
      conv_qu_pct = 100 * mean(conv_qu, na.rm = TRUE),
      mean_nll_gap_ew = mean(nll_gap_ew, na.rm = TRUE),
      mean_nll_gap_qu = mean(nll_gap_qu, na.rm = TRUE),
      .groups = "drop"
    )
  safe_write_csv(tab_conv, "Tab_Sim_Convergence_and_NLLGap.csv")
  
  invisible(list(tab_acc_1 = tab_acc_1, tab_acc_2 = tab_acc_2, tab_rt_1 = tab_rt_1, tab_rt_2 = tab_rt_2))
}

generate_sim_figures <- function(sim_df) {
  cat("Generating simulation figures...\n")
  
  acc_long <- sim_df %>%
    group_by(model, n, B) %>%
    summarise(
      rmse_ew = mean(rmse_ew, na.rm = TRUE),
      rmse_qu = mean(rmse_qu, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(c(rmse_ew, rmse_qu), names_to = "strategy", values_to = "RMSE") %>%
    mutate(
      strategy = recode(strategy, rmse_ew = "Equal Width", rmse_qu = "Equal Probability"),
      model = factor(model, levels = names(model_pretty), labels = unname(model_pretty))
    )
  
  p_rmse <- ggplot(acc_long, aes(x = B, y = RMSE, color = factor(n), linetype = strategy)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~ model, scales = "free_y") +
    scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500), labels = c("10", "20", "50", "100", "200", "500")) +
    theme_pub() +
    labs(
      x = "Number of bins (B)",
      y = "RMSE",
      color = "Sample size",
      linetype = "Binning"
    )
  save_plot(p_rmse, "Fig_Sim_RMSE_vs_B.pdf", width = 11, height = 7)
  
  rt_long <- sim_df %>%
    group_by(model, n, B) %>%
    summarise(
      time_raw = mean(time_raw, na.rm = TRUE),
      time_ew = mean(time_ew, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(B == 50) %>%
    pivot_longer(c(time_raw, time_ew), names_to = "method", values_to = "time") %>%
    mutate(
      method = recode(method, time_raw = "Raw-data MLE", time_ew = "Symbolic MLE (EW, B = 50)"),
      model = factor(model, levels = names(model_pretty), labels = unname(model_pretty))
    )
  
  p_rt <- ggplot(rt_long, aes(x = n, y = time, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ model, scales = "free_y") +
    scale_x_log10(labels = scales::comma_format()) +
    scale_y_log10() +
    theme_pub() +
    labs(
      x = "Sample size (n)",
      y = "Runtime (seconds, log scale)",
      color = "Method"
    )
  save_plot(p_rt, "Fig_Sim_Runtime_vs_n.pdf", width = 11, height = 7)
}

# ==============================
# 7. REAL-DATA INGESTION
# ==============================
find_first_matching_col <- function(df, candidates) {
  nms <- names(df)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NULL)
  nms[idx[1]]
}

load_mimic_hr <- function(csv_path) {
  df <- data.table::fread(csv_path)
  hr_col <- find_first_matching_col(df, c("heart_rate", "hr", "heartrate", "HeartRate", "value"))
  if (is.null(hr_col)) stop("Could not find a heart-rate column in: ", csv_path)
  x <- as.numeric(df[[hr_col]])
  x <- x[is.finite(x) & x > 20 & x < 250]
  if (length(x) == 0) stop("No valid heart-rate observations after filtering.")
  x
}

load_bigideas_cgm <- function(csv_path) {
  df <- data.table::fread(csv_path)
  glu_col <- find_first_matching_col(df, c("glucose", "cgm", "glucose_value", "glucose_mgdl", "value"))
  if (is.null(glu_col)) stop("Could not find a glucose column in: ", csv_path)
  x <- as.numeric(df[[glu_col]])
  x <- x[is.finite(x) & x >= 40 & x <= 400]
  if (length(x) == 0) stop("No valid glucose observations after filtering.")
  x
}

build_descriptives <- function(vec, label) {
  tibble(
    Dataset = label,
    N_Observations = length(vec),
    Min = min(vec),
    Q1 = as.numeric(quantile(vec, 0.25, na.rm = TRUE)),
    Median = median(vec, na.rm = TRUE),
    Q3 = as.numeric(quantile(vec, 0.75, na.rm = TRUE)),
    Max = max(vec, na.rm = TRUE)
  )
}

# ==============================
# 8. STREAMING / WINDOW-WISE ENGINE
# ==============================
init_hist_stream <- function(y_window, boundaries) {
  B <- length(boundaries) - 1L
  idx <- find_bin_index(y_window, boundaries)
  counts <- tabulate(idx, nbins = B)
  list(counts = as.integer(counts), boundaries = boundaries, B = B, n = length(y_window))
}

update_hist_stream <- function(hist_state, outgoing, incoming) {
  idx_out <- find_bin_index(outgoing, hist_state$boundaries)
  idx_in  <- find_bin_index(incoming, hist_state$boundaries)
  hist_state$counts[idx_out] <- hist_state$counts[idx_out] - 1L
  hist_state$counts[idx_in]  <- hist_state$counts[idx_in] + 1L
  hist_state
}

fit_symbolic_from_hist <- function(hist_obj, model, init_par = NULL) {
  if (is.null(init_par)) {
    mids <- 0.5 * (hist_obj$boundaries[-1] + hist_obj$boundaries[-length(hist_obj$boundaries)])
    pseudo_y <- rep(mids, hist_obj$counts)
    init_par <- get_init(pseudo_y, model)
  }
  obj <- function(z) symbolic_nll(hist_obj, model, z)
  t0 <- proc.time()["elapsed"]
  opt <- tryCatch(
    suppressWarnings(optim(init_par, obj, method = "BFGS", control = list(maxit = OPT_MAXIT))),
    error = function(e) list(par = init_par, value = NA_real_, convergence = 99, message = conditionMessage(e))
  )
  elapsed <- proc.time()["elapsed"] - t0
  list(
    est = transform_params(opt$par, model),
    nll = opt$value,
    converged = isTRUE(opt$convergence == 0),
    time = as.numeric(elapsed),
    par_unconstrained = opt$par
  )
}

candidate_models_for_dataset <- function(dataset_name) {
  if (dataset_name == "MIMIC-III") return(c("normal", "gamma"))
  if (dataset_name == "BIG IDEAs") return(c("gamma", "lognormal"))
  stop("Unknown dataset name.")
}

run_real_data_windows <- function(stream, dataset_name, win_size, step, B, cfg) {
  models <- candidate_models_for_dataset(dataset_name)
  ends <- seq(win_size, length(stream), by = step)
  if (is.finite(cfg$stream_max_win)) ends <- head(ends, cfg$stream_max_win)
  n_windows <- length(ends)
  if (n_windows < 1) stop("No windows available for dataset: ", dataset_name)
  
  cat(sprintf("Running %s windows (B = %d)...\n", dataset_name, B))
  
  boundaries <- seq(min(stream), max(stream), length.out = B + 1)
  boundaries[1] <- boundaries[1] - 1e-8
  boundaries[length(boundaries)] <- boundaries[length(boundaries)] + 1e-8
  
  initial_window <- stream[1:win_size]
  hist_state <- init_hist_stream(initial_window, boundaries)
  
  results <- vector("list", n_windows * length(models))
  counter <- 1L
  prev_raw_par <- list()
  prev_sym_par <- list()
  
  bar_width <- 50
  
  for (w_idx in seq_along(ends)) {
    end_idx <- ends[w_idx]
    start_idx <- end_idx - win_size + 1L
    y_win <- stream[start_idx:end_idx]
    
    if (w_idx == 1L) {
      hist_state <- init_hist_stream(y_win, boundaries)
    } else {
      prev_end <- ends[w_idx - 1L]
      new_points <- stream[(prev_end + 1L):end_idx]
      old_start <- prev_end - win_size + 1L
      old_points <- stream[old_start:(old_start + length(new_points) - 1L)]
      for (k in seq_along(new_points)) {
        hist_state <- update_hist_stream(hist_state, outgoing = old_points[k], incoming = new_points[k])
      }
    }
    
    hist_obj <- list(boundaries = hist_state$boundaries, counts = hist_state$counts, B = hist_state$B, n = hist_state$n)
    
    for (model in models) {
      raw_fit <- fit_mle(y_win, model, init = prev_raw_par[[model]])
      sym_fit <- fit_symbolic_from_hist(hist_obj, model, init_par = prev_sym_par[[model]])
      
      # Fixed unname() wrapper to prevent shape.shape/mu.mu recursion bugs
      prev_raw_par[[model]] <- switch(
        model,
        normal    = unname(c(raw_fit$est["mu"], log(raw_fit$est["sigma"]))),
        exp       = unname(c(log(raw_fit$est["rate"]))),
        gamma     = unname(c(log(raw_fit$est["shape"]), log(raw_fit$est["rate"]))),
        weibull   = unname(c(log(raw_fit$est["shape"]), log(raw_fit$est["scale"]))),
        lognormal = unname(c(raw_fit$est["meanlog"], log(raw_fit$est["sdlog"])))
      )
      prev_sym_par[[model]] <- unname(sym_fit$par_unconstrained)
      
      param_names <- union(names(raw_fit$est), names(sym_fit$est))
      p1 <- param_names[1]
      p2 <- if (length(param_names) >= 2) param_names[2] else NA_character_
      
      results[[counter]] <- tibble(
        Dataset = dataset_name,
        Window = w_idx,
        WindowStart = start_idx,
        WindowEnd = end_idx,
        Model = model,
        Raw_Time = raw_fit$time,
        Sym_Time = sym_fit$time,
        Raw_NLL = raw_fit$nll,
        Sym_NLL = sym_fit$nll,
        SpeedUp = safe_divide(raw_fit$time, sym_fit$time),
        Param1_Name = p1,
        Param1_Raw = unname(raw_fit$est[p1]),
        Param1_Sym = unname(sym_fit$est[p1]),
        Param2_Name = p2,
        Param2_Raw = if (!is.na(p2)) unname(raw_fit$est[p2]) else NA_real_,
        Param2_Sym = if (!is.na(p2)) unname(sym_fit$est[p2]) else NA_real_,
        Param_AbsDiff = param_absdiff(raw_fit$est, sym_fit$est),
        Raw_Conv = raw_fit$converged,
        Sym_Conv = sym_fit$converged
      )
      counter <- counter + 1L
    }
    
    # Progress bar monitor
    pct <- w_idx / n_windows
    filled <- round(pct * bar_width)
    empty <- bar_width - filled
    bar <- paste0(rep("=", filled), collapse = "")
    if (empty > 0) {
      bar <- paste0(bar, ">")
      empty <- empty - 1
    }
    space <- paste0(rep("-", max(0, empty)), collapse = "")
    
    cat(sprintf("\r  Progress: [%s%s] %d / %d (%5.1f%%)", bar, space, w_idx, n_windows, 100 * pct))
    flush.console()
  }
  cat("\n") # New line after 100% completion
  
  bind_rows(results)
}

select_best_model <- function(window_df) {
  window_df %>%
    group_by(Model) %>%
    summarise(mean_raw_nll = mean(Raw_NLL, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_raw_nll) %>%
    slice(1) %>%
    pull(Model)
}

summarise_realdata_runtime <- function(df) {
  best_model <- select_best_model(df)
  best_df <- df %>% filter(Model == best_model)
  tibble(
    preferred_model = best_model,
    raw_median = median(best_df$Raw_Time, na.rm = TRUE),
    raw_iqr = IQR(best_df$Raw_Time, na.rm = TRUE),
    sym_median = median(best_df$Sym_Time, na.rm = TRUE),
    sym_iqr = IQR(best_df$Sym_Time, na.rm = TRUE),
    speedup_median = median(best_df$SpeedUp, na.rm = TRUE),
    speedup_iqr = IQR(best_df$SpeedUp, na.rm = TRUE),
    absdiff_median = median(best_df$Param_AbsDiff, na.rm = TRUE),
    absdiff_max = max(best_df$Param_AbsDiff, na.rm = TRUE),
    n_windows = nrow(best_df)
  )
}

run_sensitivity_grid <- function(stream, dataset_name, B_grid, W_grid, step, cfg) {
  out <- list()
  idx <- 1L
  for (B in B_grid) {
    for (W in W_grid) {
      tmp <- run_real_data_windows(stream, dataset_name, win_size = W, step = step, B = B, cfg = cfg)
      best_model <- select_best_model(tmp)
      tmp_best <- tmp %>% filter(Model == best_model)
      out[[idx]] <- tibble(
        Dataset = dataset_name,
        B = B,
        WindowSize = W,
        PreferredModel = best_model,
        Avg_Sym_Runtime = mean(tmp_best$Sym_Time, na.rm = TRUE),
        Median_SpeedUp = median(tmp_best$SpeedUp, na.rm = TRUE),
        Discrepancy = median(tmp_best$Param_AbsDiff, na.rm = TRUE)
      )
      idx <- idx + 1L
    }
  }
  bind_rows(out)
}

# ==============================
# 9. REAL-DATA OUTPUTS
# ==============================
generate_realdata_outputs <- function(cfg) {
  cat("Running real-data validation...\n")
  
  hr_file  <- file.path(mainDir, "data", "mimic_hr.csv")
  cgm_file <- file.path(mainDir, "data", "bigideas_cgm.csv")
  
  if (!file.exists(hr_file)) stop("Missing required file: ", hr_file)
  if (!file.exists(cgm_file)) stop("Missing required file: ", cgm_file)
  
  hr_data  <- load_mimic_hr(hr_file)
  cgm_data <- load_bigideas_cgm(cgm_file)
  
  tab_mimic_desc <- build_descriptives(hr_data, "MIMIC-III Heart Rate Stream")
  tab_bigideas_desc <- build_descriptives(cgm_data, "BIG IDEAs CGM Glucose Stream")
  
  safe_write_csv(tab_mimic_desc, "Tab_MIMIC_Desc.csv")
  safe_write_csv(tab_bigideas_desc, "Tab_BIGIDEAs_Desc.csv")
  
  mimic_res <- run_real_data_windows(
    stream = hr_data,
    dataset_name = "MIMIC-III",
    win_size = cfg$stream_win_hr,
    step = cfg$stream_step_hr,
    B = cfg$realdata_B,
    cfg = cfg
  )
  
  bigideas_res <- run_real_data_windows(
    stream = cgm_data,
    dataset_name = "BIG IDEAs",
    win_size = cfg$stream_win_cgm,
    step = cfg$stream_step_cgm,
    B = cfg$realdata_B,
    cfg = cfg
  )
  
  safe_write_csv(mimic_res, "realdata_mimic_window_results.csv", folder = "results/realdata")
  safe_write_csv(bigideas_res, "realdata_bigideas_window_results.csv", folder = "results/realdata")
  
  mimic_sum <- summarise_realdata_runtime(mimic_res)
  bigideas_sum <- summarise_realdata_runtime(bigideas_res)
  
  tab_real_runtime <- tibble(
    Dataset = c("MIMIC-III", "BIG IDEAs"),
    Preferred_Model = c(mimic_sum$preferred_model, bigideas_sum$preferred_model),
    `Raw MLE (s)` = c(mimic_sum$raw_median, bigideas_sum$raw_median),
    `Symbolic MLE (s)` = c(mimic_sum$sym_median, bigideas_sum$sym_median),
    `Speed-up` = c(mimic_sum$speedup_median, bigideas_sum$speedup_median),
    `Speed-up IQR` = c(mimic_sum$speedup_iqr, bigideas_sum$speedup_iqr),
    `Median Abs Param Diff` = c(mimic_sum$absdiff_median, bigideas_sum$absdiff_median),
    `Max Abs Param Diff` = c(mimic_sum$absdiff_max, bigideas_sum$absdiff_max),
    `N Windows` = c(mimic_sum$n_windows, bigideas_sum$n_windows)
  )
  safe_write_csv(tab_real_runtime, "Tab_RealData_Runtime.csv")
  
  mimic_best <- mimic_res %>% filter(Model == mimic_sum$preferred_model)
  bigideas_best <- bigideas_res %>% filter(Model == bigideas_sum$preferred_model)
  
  mimic_plot_df <- bind_rows(
    mimic_best %>%
      transmute(Window, Panel = paste0("Parameter: ", Param1_Name), Method = "Raw", Value = Param1_Raw),
    mimic_best %>%
      transmute(Window, Panel = paste0("Parameter: ", Param1_Name), Method = "Symbolic", Value = Param1_Sym),
    mimic_best %>%
      transmute(Window, Panel = "Per-window runtime", Method = "Raw", Value = Raw_Time),
    mimic_best %>%
      transmute(Window, Panel = "Per-window runtime", Method = "Symbolic", Value = Sym_Time)
  )
  
  p_mimic <- ggplot(mimic_plot_df, aes(x = Window, y = Value, color = Method, linetype = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
    theme_pub() +
    labs(
      x = "Window index",
      y = NULL,
      color = "Estimator",
      linetype = "Estimator"
    )
  save_plot(p_mimic, "Fig_MIMIC_Window_Params.pdf", width = 10, height = 7)
  
  bigideas_plot_df <- bind_rows(
    bigideas_best %>%
      transmute(Window, Panel = paste0("Parameter: ", Param1_Name), Method = "Raw", Value = Param1_Raw),
    bigideas_best %>%
      transmute(Window, Panel = paste0("Parameter: ", Param1_Name), Method = "Symbolic", Value = Param1_Sym),
    bigideas_best %>%
      transmute(Window, Panel = "Per-window runtime", Method = "Raw", Value = Raw_Time),
    bigideas_best %>%
      transmute(Window, Panel = "Per-window runtime", Method = "Symbolic", Value = Sym_Time)
  )
  
  p_bigideas <- ggplot(bigideas_plot_df, aes(x = Window, y = Value, color = Method, linetype = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
    theme_pub() +
    labs(
      x = "Window index",
      y = NULL,
      color = "Estimator",
      linetype = "Estimator"
    )
  save_plot(p_bigideas, "Fig_BIGIDEAs_Window_Params.pdf", width = 10, height = 7)
  
  sens_mimic <- run_sensitivity_grid(
    stream = hr_data,
    dataset_name = "MIMIC-III",
    B_grid = cfg$sens_B_grid,
    W_grid = cfg$sens_W_grid_hr,
    step = cfg$stream_step_hr,
    cfg = cfg
  )
  sens_bigideas <- run_sensitivity_grid(
    stream = cgm_data,
    dataset_name = "BIG IDEAs",
    B_grid = cfg$sens_B_grid,
    W_grid = cfg$sens_W_grid_cgm,
    step = cfg$stream_step_cgm,
    cfg = cfg
  )
  
  sens_df <- bind_rows(sens_mimic, sens_bigideas)
  safe_write_csv(sens_df, "realdata_sensitivity_grid.csv", folder = "results/realdata")
  
  sens_plot_df <- sens_df %>%
    pivot_longer(cols = c(Avg_Sym_Runtime, Discrepancy), names_to = "Metric", values_to = "Value") %>%
    mutate(
      Metric = recode(
        Metric,
        Avg_Sym_Runtime = "Average symbolic runtime (seconds)",
        Discrepancy = "Median absolute parameter discrepancy"
      )
    )
  
  p_sens <- ggplot(sens_plot_df, aes(x = B, y = Value, color = factor(WindowSize), group = WindowSize)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_grid(Metric ~ Dataset, scales = "free_y") +
    theme_pub() +
    labs(
      x = "Number of bins (B)",
      y = NULL,
      color = "Window size"
    )
  save_plot(p_sens, "Fig_RealData_Sensitivity.pdf", width = 11, height = 7)
  
  invisible(list(
    mimic_res = mimic_res,
    bigideas_res = bigideas_res,
    sens_df = sens_df,
    tab_real_runtime = tab_real_runtime
  ))
}

# ==============================
# 10. MANUSCRIPT-READY EXPORT HELPERS
# ==============================
make_manuscript_summary <- function(sim_df, real_out) {
  tab_acc <- sim_df %>%
    group_by(model, n, B) %>%
    summarise(
      rel_rmse_ew = mean(rel_rmse_ew, na.rm = TRUE),
      rel_rmse_qu = mean(rel_rmse_qu, na.rm = TRUE),
      speedup_ew = mean(speedup_ew, na.rm = TRUE),
      speedup_qu = mean(speedup_qu, na.rm = TRUE),
      .groups = "drop"
    )
  
  largest_n <- max(sim_df$n)
  speed_rng <- tab_acc %>%
    filter(n == largest_n) %>%
    summarise(min_speed = min(speedup_ew, na.rm = TRUE), max_speed = max(speedup_ew, na.rm = TRUE))
  
  rel_rmse_max_B100 <- tab_acc %>%
    filter(n == largest_n, B >= 100) %>%
    summarise(max_rel = max(rel_rmse_ew, na.rm = TRUE))
  
  summary_text <- tibble(
    item = c(
      "Largest-n equal-width speed-up range",
      "Largest-n max relative RMSE for EW when B >= 100",
      "MIMIC preferred model",
      "MIMIC median speed-up",
      "BIG IDEAs preferred model",
      "BIG IDEAs median speed-up"
    ),
    value = c(
      sprintf("%.2f-fold to %.2f-fold", speed_rng$min_speed, speed_rng$max_speed),
      sprintf("%.3f", rel_rmse_max_B100$max_rel),
      real_out$tab_real_runtime$Preferred_Model[real_out$tab_real_runtime$Dataset == "MIMIC-III"],
      sprintf("%.3f", real_out$tab_real_runtime$`Speed-up`[real_out$tab_real_runtime$Dataset == "MIMIC-III"]),
      real_out$tab_real_runtime$Preferred_Model[real_out$tab_real_runtime$Dataset == "BIG IDEAs"],
      sprintf("%.3f", real_out$tab_real_runtime$`Speed-up`[real_out$tab_real_runtime$Dataset == "BIG IDEAs"])
    )
  )
  
  safe_write_csv(summary_text, "Tab_Manuscript_Summary_Stats.csv")
}

# ==============================
# 11. MASTER CONTROLLER
# ==============================
run_full_pipeline <- function(cfg = ACTIVE_CONFIG) {
  pipeline_start_time <- Sys.time() # START TIMER
  
  sim_df <- run_simulation_study(cfg)
  generate_sim_tables(sim_df)
  generate_sim_figures(sim_df)
  real_out <- generate_realdata_outputs(cfg)
  make_manuscript_summary(sim_df, real_out)
  
  pipeline_end_time <- Sys.time() # END TIMER
  elapsed_secs <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "secs"))
  
  hrs <- floor(elapsed_secs / 3600)
  mins <- floor((elapsed_secs %% 3600) / 60)
  secs <- round(elapsed_secs %% 60, 1)
  
  time_str <- if (hrs > 0) {
    sprintf("%d hours, %d minutes, %.1f seconds", hrs, mins, secs)
  } else if (mins > 0) {
    sprintf("%d minutes, %.1f seconds", mins, secs)
  } else {
    sprintf("%.1f seconds", secs)
  }
  
  cat("\n====================================================================\n")
  cat(" FULL PIPELINE COMPLETE\n")
  cat(sprintf(" Total Execution Time: %s\n", time_str))
  cat(" Tables saved to: ", file.path(mainDir, "tables"), "\n", sep = "")
  cat(" Figures saved to: ", file.path(mainDir, "figures"), "\n", sep = "")
  cat(" Results saved to: ", file.path(mainDir, "results"), "\n", sep = "")
  cat("====================================================================\n")
  
  invisible(list(simulation = sim_df, realdata = real_out))
}

# ==============================
# 12. EXECUTE
# ==============================
run_full_pipeline(ACTIVE_CONFIG)

