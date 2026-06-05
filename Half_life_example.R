# ==========================================================================
# AUTOCORRELATION-BASED HYDROLOGICAL MEMORY METRICS
# ==========================================================================
#
# Purpose:
#   Quantify catchment-scale hydrological memory from daily streamflow
#   time series using autocorrelation functions (ACFs).
#
# Workflow:
#   1. Generate or load daily streamflow time series.
#   2. Calculate autocorrelation functions across a range of lag times.
#   3. Fit an exponential decay model to the positive portion of the ACF:
#
#          r(lag) = exp(-lag / tau)
#
#   4. Estimate model-based memory metrics:
#          tau        = characteristic decay timescale
#          half_life  = lag where fitted autocorrelation reaches 0.5
#          r2         = goodness-of-fit of exponential decay model
#
#   5. Estimate empirical memory metrics from the raw ACF:
#          lag_cross0 = interpolated lag where ACF reaches zero
#          lag_at_eps = interpolated lag where ACF falls below 0.05
#
# Outputs:
#   result_long   = autocorrelation function for each catchment and lag
#   final_metrics = summary memory metrics for each catchment
#
# Interpretation:
#   Larger values of tau, half_life, lag_cross0, and lag_at_eps indicate
#   greater persistence and longer hydrological memory.
#
# ==========================================================================
library(data.table)
library(ggplot2)

# ==========================================================================
# 1. CREATE SYNTHETIC DAILY STREAMFLOW DATA
# ==========================================================================
#
# Inputs:
#   n_years      = length of simulation (years)
#   n_days       = total number of daily observations
#   n_ids        = number of synthetic catchments
#
# Outputs:
#   test         = data.table containing:
#                    date     = daily timestamps
#                    ID_*     = synthetic daily streamflow series
#
# Purpose:
#   Generate synthetic streamflow records with varying persistence and
#   seasonality for testing the autocorrelation-based memory metrics.
#
# ==========================================================================

set.seed(123)

n_years <- 10
n_days <- n_years * 365
n_ids <- 20

dates <- seq.Date(
  from = as.Date("2000-01-01"),
  by = "day",
  length.out = n_days
)

test <- data.table(date = dates)

for (i in seq_len(n_ids)) {
  
  ID <- paste0("ID_", sprintf("%03d", i))
  
  # Persistence parameter: higher phi = longer hydrological memory.
  phi <- runif(1, 0.45, 0.98)
  
  # Random noise.
  eps <- rnorm(n_days, mean = 0, sd = 1)
  
  # AR(1)-type storage/memory process.
  x <- numeric(n_days)
  x[1] <- eps[1]
  
  for (t in 2:n_days) {
    x[t] <- phi * x[t - 1] + eps[t]
  }
  
  # Seasonal component.
  seasonal <- sin(2 * pi * seq_len(n_days) / 365)
  
  # Positive streamflow-like values.
  Q <- exp(1 + 0.15 * x + 0.25 * seasonal)
  
  test[, (ID) := Q]
}

# Inspect synthetic streamflow table.
test[]


# ==========================================================================
# 2. COMPUTE AUTOCORRELATION FUNCTIONS
# ==========================================================================
#
# Inputs:
#   test         = daily streamflow table
#   lags         = vector of lag times (days)
#
# Outputs:
#   result_long  = long-format autocorrelation table containing:
#
#                  ID         = catchment identifier
#                  lag_days   = separation between paired observations used to compute
#            autocorrelation, e.g. lag_days = 30 corresponds to
#            Corr(Q_t, Q_(t+30)) calculated across the full record.
#                  r          = autocorrelation at lag_days
#                  n_pairs    = number of valid lagged pairs
#                  ci_lower   = lower 95% confidence interval
#                  ci_upper   = upper 95% confidence interval
#
# Purpose:
# Quantify the average dependence between streamflow at time t and
# streamflow at time t + lag across the full record.
#
# ==========================================================================

# Daily lags.
lags <- seq(1, 365, 10)

# ID columns only, excluding date.
id_cols <- names(test)[-1]

# Convert streamflow columns to matrix.
X <- as.matrix(test[, ..id_cols])


# Correlation and Fisher-z confidence interval at one lag for one series.
acf_at_lag_with_ci <- function(x, L, conf.level = 0.95) {
  
  n <- length(x)
  
  if (n <= L) {
    return(list(r = NA_real_, n = 0L, lo = NA_real_, hi = NA_real_))
  }
  
  x1 <- x[(L + 1):n]
  x0 <- x[1:(n - L)]
  
  ok <- is.finite(x1) & is.finite(x0)
  x1 <- x1[ok]
  x0 <- x0[ok]
  
  m <- length(x1)
  
  if (m < 4L) {
    return(list(r = NA_real_, n = m, lo = NA_real_, hi = NA_real_))
  }
  
  r <- suppressWarnings(cor(x1, x0))
  
  if (!is.finite(r)) {
    return(list(r = NA_real_, n = m, lo = NA_real_, hi = NA_real_))
  }
  
  # Fisher z CI.
  r_clamped <- max(min(r, 0.999999), -0.999999)
  z <- atanh(r_clamped)
  se <- 1 / sqrt(m - 3)
  zcrit <- qnorm(1 - (1 - conf.level) / 2)
  
  lo <- tanh(z - zcrit * se)
  hi <- tanh(z + zcrit * se)
  
  list(r = r, n = m, lo = lo, hi = hi)
}


# Compute for all IDs and lags.
res_list <- lapply(seq_len(ncol(X)), function(j) {
  
  id <- colnames(X)[j]
  v <- X[, j]
  
  lapply(lags, function(L) {
    
    out <- acf_at_lag_with_ci(v, L)
    
    data.table(
      ID       = id,
      lag_days = L,
      r        = out$r,
      n_pairs  = out$n,
      ci_lower = out$lo,
      ci_upper = out$hi
    )
  })
})

result_long <- rbindlist(unlist(res_list, recursive = FALSE))

# This is equivalent to your precomputed auto_correlation.csv table.
dt_all <- copy(result_long)

dt_all[]


# ==========================================================================
# 3. ESTIMATE CHARACTERISTIC DECAY TIME (TAU)
# ==========================================================================
#
# Inputs:
#   dt_all       = autocorrelation table produced above
#
# Outputs:
#   decay_est    = memory metrics derived from exponential decay fitting:
#
#                  tau        = characteristic decay time (days)
#                  tau_lo     = lower 95% confidence interval for tau
#                  tau_hi     = upper 95% confidence interval for tau
#                  half_life  = lag where fitted ACF reaches 50%
#                  slope      = fitted log-linear decay slope
#                  # r2 = goodness-of-fit for the through-origin log-linear model.
#      Calculated using SST = Σwy² because the regression is forced
#      through the origin (log(r)=0 at lag=0).
#                  n_pts      = number of positive ACF points used
#
# Method:
#   Fit:
#
#       r(lag) = exp(-lag / tau)
#
#   by estimating a weighted through-origin regression:
#
#       log(r) = slope × lag
#
#   using only positive autocorrelation values.
#
# Purpose:
#   Summarize catchment memory as a single characteristic decay time.
# Note: Check the r2 before trusting tau. 
#
# ==========================================================================

dt_all[, `:=`(
  lag_days = as.numeric(lag_days),
  r        = as.numeric(r),
  n_pairs  = as.numeric(n_pairs)
)]

setorder(dt_all, ID, lag_days)

# Exponential decay fitting requires log(r).
# Therefore only positive autocorrelation values are used.
# Zero-crossing metrics are estimated separately using the full ACF.

dt_tau <- dt_all[
  is.finite(r) &
    r > 0 &
    is.finite(lag_days) &
    lag_days > 0
]

eps1 <- 1 - 1e-6

dt_tau[, log_r := log(pmin(r, eps1))]

decay_est <- dt_tau[, {
  
  x <- lag_days
  y <- log_r
  w <- fifelse(is.finite(n_pairs) & n_pairs > 0, n_pairs, 1)
  
  keep <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[keep]
  y <- y[keep]
  w <- w[keep]
  
  n_pts <- length(x)
  
  if (n_pts < 3L) {
    
    .(
      tau       = NA_real_,
      tau_lo    = NA_real_,
      tau_hi    = NA_real_,
      half_life = NA_real_,
      slope     = NA_real_,
      r2        = NA_real_,
      n_pts     = n_pts
    )
    
  } else {
    
    Sxx <- sum(w * x * x)
    Sxy <- sum(w * x * y)
    s   <- Sxy / Sxx
    
    yhat <- s * x
    sse  <- sum(w * (y - yhat)^2)
    sst  <- sum(w * y^2)
    r2   <- if (sst > 0) 1 - sse / sst else NA_real_
    
    if (!is.finite(s) || s >= 0) {
      
      .(
        tau       = NA_real_,
        tau_lo    = NA_real_,
        tau_hi    = NA_real_,
        half_life = NA_real_,
        slope     = s,
        r2        = r2,
        n_pts     = n_pts
      )
      
    } else {
      
      tau <- -1 / s
      half_life <- tau * log(2)
      
      sigma2 <- if (n_pts > 1) sse / (n_pts - 1) else NA_real_
      se_s   <- sqrt(sigma2 / Sxx)
      se_tau <- se_s / (s^2)
      z      <- qnorm(0.975)
      
      .(
        tau       = tau,
        tau_lo    = tau - z * se_tau,
        tau_hi    = tau + z * se_tau,
        half_life = half_life,
        slope     = s,
        r2        = r2,
        n_pts     = n_pts
      )
    }
  }
  
}, by = ID]


# ==========================================================================
# 4. ESTIMATE ZERO-CROSSING AND NEAR-ZERO MEMORY LAGS
# ==========================================================================
#
# Inputs:
#   dt_all       = full autocorrelation table
#   eps0         = practical near-zero threshold (default = 0.05)
#
# Outputs:
#   zero_cross   = empirical memory metrics:
#
#                  lag_cross0      = interpolated zero crossing (days)
#                  lag_first_le0   = first observed lag with r <= 0
#                  lag_at_eps      = interpolated lag where r <= eps0
#                  min_r           = minimum observed autocorrelation
#
# Method:
#   Uses the raw autocorrelation function (including negative values)
#   and linearly interpolates between neighboring lags.
#
# Purpose:
#   Estimate the point where autocorrelation effectively disappears,
#   without assuming an exponential decay model.
#
# ==========================================================================

eps0 <- 0.05

zero_cross <- dt_all[, {
  
  x <- lag_days
  y <- r
  
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  idx_le0 <- which(y <= 0)
  lag_first_le0 <- if (length(idx_le0)) x[idx_le0[1]] else NA_real_
  
  lag_cross0 <- NA_real_
  
  if (length(idx_le0)) {
    
    i2 <- idx_le0[1]
    
    if (i2 == 1L || !is.finite(y[i2 - 1])) {
      
      lag_cross0 <- x[i2]
      
    } else {
      
      i1 <- i2 - 1L
      y1 <- y[i1]
      y2 <- y[i2]
      
      lag_cross0 <- if (
        is.finite(y1) &&
        is.finite(y2) &&
        (y2 - y1) != 0
      ) {
        x[i1] + (0 - y1) * (x[i2] - x[i1]) / (y2 - y1)
      } else {
        x[i2]
      }
    }
  }
  
  lag_at_eps <- {
    
    idx_eps <- which(y <= eps0)
    
    if (!length(idx_eps)) {
      
      NA_real_
      
    } else {
      
      i2 <- idx_eps[1]
      
      if (i2 == 1L) {
        
        x[i2]
        
      } else {
        
        i1 <- i2 - 1L
        y1 <- y[i1]
        y2 <- y[i2]
        
        if (
          is.finite(y1) &&
          is.finite(y2) &&
          (y2 - y1) != 0
        ) {
          x[i1] + (eps0 - y1) * (x[i2] - x[i1]) / (y2 - y1)
        } else {
          x[i2]
        }
      }
    }
  }
  
  .(
    lag_cross0    = lag_cross0,
    lag_first_le0 = lag_first_le0,
    lag_at_eps    = lag_at_eps,
    min_r         = if (length(y)) min(y, na.rm = TRUE) else NA_real_
  )
  
}, by = ID]


# ==========================================================================
# 5. MERGE MEMORY METRICS
# ==========================================================================
#
# Inputs:
#   decay_est    = model-based memory metrics
#   zero_cross   = empirical memory metrics
#
# Outputs:
#   final_metrics = one row per catchment containing:
#
#                   tau
#                   tau_lo
#                   tau_hi
#                   half_life
#                   slope
#                   r2
#                   n_pts
#                   lag_cross0
#                   lag_first_le0
#                   lag_at_eps
#                   min_r
#
# Purpose:
#   Combine all memory metrics into a single summary table for
#   subsequent analysis, visualization, or comparison between
#   catchments and hydrological regimes.
#
# ==========================================================================

final_metrics <- zero_cross[decay_est, on = "ID"]

final_metrics[]


# ==========================================================================
# INTERPRETATION OF MEMORY METRICS
# ==========================================================================
#
# tau
#   Characteristic exponential decay timescale.
#   At lag = tau, the fitted autocorrelation has declined to exp(-1)
#   ≈ 0.368 of its initial value.
#
# half_life
#   Lag where fitted autocorrelation declines to 50%.
#
# Relationship:
#   half_life = tau * log(2)
#
# Therefore both metrics describe the same fitted exponential decay
# curve using different reference levels.
#
# lag_at_eps
#   Lag where autocorrelation first falls below 0.05.
#   Represents practical loss of memory.
#
# lag_cross0
#   Lag where autocorrelation becomes zero.
#
# r2
#   Indicates how closely the ACF follows exponential decay.
#
# Example:
#
#   Catchment A:
#       tau = 10 days
#       half_life = 7 days
#       lag_at_eps = 25 days
#
#   Catchment B:
#       tau = 60 days
#       half_life = 42 days
#       lag_at_eps = 120 days
#
#   Catchment B exhibits substantially longer hydrological memory.
#
# ==========================================================================

# ==========================================================================
# 6. OPTIONAL PLOTS
# ==========================================================================

ggplot(result_long, aes(x = lag_days, y = r, group = ID)) +
  geom_line(alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = eps0, linetype = "dotted") +
  theme_bw() +
  labs(
    x = "Lag in days",
    y = "Autocorrelation",
    title = "Synthetic autocorrelation functions"
  )


