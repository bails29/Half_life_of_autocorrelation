# Hydrological Memory Metrics from Streamflow Autocorrelation

## Overview

This workflow quantifies catchment-scale hydrological memory using autocorrelation functions (ACFs) derived from daily streamflow time series.

The approach computes empirical autocorrelation curves over a range of lag times and summarizes these curves using both model-based and empirical memory metrics. The resulting metrics provide complementary measures of how long streamflow anomalies persist within a catchment.

The workflow was designed for comparing hydrological memory across large samples of catchments.

---

## Method Summary

### Step 1: Compute autocorrelation functions

For each catchment, autocorrelation is calculated for a range of lag times:

[
r(L) = \mathrm{Corr}(Q_t, Q_{t+L})
]

where:

* (Q_t) is streamflow at time (t)
* (L) is lag time (days)
* (r(L)) is the autocorrelation at lag (L)

Autocorrelation is calculated using all valid observation pairs separated by a given lag.

For example, a lag of 30 days measures the correlation between flow on day (t) and flow on day (t+30) across the entire record.

The output of this step is a full autocorrelation function (ACF) for each catchment.

---

### Step 2: Fit an exponential decay model

The positive portion of each ACF is approximated using:

[
r(L)=e^{-L/\tau}
]

where:

* (r(L)) = autocorrelation
* (L) = lag (days)
* (\tau) = characteristic decay timescale

Taking logarithms gives:

[
\log(r)=sL
]

where:

[
\tau=-\frac{1}{s}
]

A weighted through-origin regression is fitted using the number of valid observation pairs at each lag as weights.

Only positive autocorrelation values are used because (\log(r)) is undefined for (r \le 0).

---

### Step 3: Calculate empirical memory metrics

Additional metrics are estimated directly from the observed ACF without assuming exponential decay.

These metrics are based on interpolation between neighboring ACF values and provide alternative measures of memory duration.

---

## Input Data

The workflow expects a table containing:

| Column | Description                                |
| ------ | ------------------------------------------ |
| date   | Daily timestamps                           |
| ID_*   | Daily streamflow values for each catchment |

Example:

| date       | ID_001 | ID_002 |
| ---------- | ------ | ------ |
| 2000-01-01 | 12.4   | 4.7    |
| 2000-01-02 | 11.8   | 5.1    |
| ...        | ...    | ...    |

Each catchment is stored in a separate column.

---

## Outputs

### 1. Autocorrelation table (`result_long`)

Contains the complete ACF for every catchment.

| Variable | Description                   |
| -------- | ----------------------------- |
| ID       | Catchment identifier          |
| lag_days | Lag length (days)             |
| r        | Autocorrelation               |
| n_pairs  | Number of paired observations |
| ci_lower | Lower 95% confidence interval |
| ci_upper | Upper 95% confidence interval |

Each row corresponds to one catchment and one lag.

---

### 2. Memory metrics table (`final_metrics`)

Contains one row per catchment.

#### Model-based metrics

##### tau

Characteristic exponential decay timescale:

[
r(\tau)=e^{-1}\approx0.368
]

Larger values indicate slower decay of autocorrelation and greater hydrological memory.

---

##### half_life

Lag at which the fitted exponential decay reaches 50%:

[
\mathrm{half_life}=\tau\ln(2)
]

Because half-life and tau are derived from the same fitted model, they contain equivalent information expressed on different scales.

---

##### r2

Goodness-of-fit of the exponential decay model.

Higher values indicate that the ACF closely follows a simple exponential decay.

Low values may indicate multiple characteristic timescales, strong seasonality, or non-exponential behavior.

---

##### slope

Slope of the fitted relationship:

[
\log(r)=sL
]

More negative slopes correspond to faster memory loss.

---

##### tau_lo, tau_hi

Approximate 95% confidence interval for tau estimated using the delta method.

---

#### Empirical metrics

##### lag_cross0

Interpolated lag at which the observed ACF reaches zero.

Represents the point where positive persistence disappears.

---

##### lag_first_le0

First sampled lag where autocorrelation is less than or equal to zero.

Unlike `lag_cross0`, this value is not interpolated.

---

##### lag_at_eps

Interpolated lag where autocorrelation falls below a threshold value of:

[
r = 0.05
]

This can be interpreted as a practical estimate of memory duration.

---

##### min_r

Minimum observed autocorrelation value across all evaluated lags.

---

## Interpretation

All memory metrics increase with hydrological persistence.

For example:

| Catchment | tau     | half_life | lag_at_eps |
| --------- | ------- | --------- | ---------- |
| A         | 10 days | 7 days    | 25 days    |
| B         | 60 days | 42 days   | 120 days   |

Catchment B exhibits substantially longer hydrological memory than Catchment A.

In general:

* Large tau → slow exponential decay
* Large half_life → long persistence
* Large lag_at_eps → long practical memory
* Large lag_cross0 → autocorrelation remains positive for longer

---

## Notes

* Autocorrelation is calculated across the full record and is not conditioned on season, flow state, or event type.
* The resulting metrics represent average catchment-scale memory.
* Exponential decay fitting uses only positive autocorrelation values.
* Empirical metrics are calculated from the full ACF, including negative values.
* The goodness-of-fit statistic (`r2`) should be examined before interpreting tau as a meaningful single-timescale description of catchment memory.
