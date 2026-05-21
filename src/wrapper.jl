"""
    fit_single_event(left, right, event, zzi, timepoints; niter=1000, tol=1e-4, factor=0.75)

Fit a proportional hazards model for interval-censored single-event data.

Estimates regression coefficients in a proportional hazards model for
interval-censored time-to-event data. Each individual is observed within an
interval `(left, right]`, with an indicator specifying whether the event
occurred. Computation is performed by the multistate estimation routine
`fit_multistate`.

# Arguments

- `left::Vector{Float64}`: Left endpoints of observation intervals. For the
  `i`-th individual, `left[i]` is the lower bound of the interval in which the
  event time lies. In cases of right censoring, `left[i]` is the last follow-up
  time.
- `right::Vector{Float64}`: Right endpoints of observation intervals. For the
  `i`-th individual, `right[i]` is the upper bound of the interval in which the
  event time lies.
- `event::Vector{Bool}`: Event status. `event[i] == true` indicates that the
  event occurred within `(left[i], right[i]]`; `event[i] == false` indicates
  right censoring.
- `zzi::Matrix{Float64}`: Covariate matrix. The `i`-th row corresponds to the
  `i`-th individual and each column to a predictor variable.
- `timepoints::Vector{Float64}`: Timepoints for the discrete approximation.
  Values must be finite, positive, and strictly increasing. `maximum(timepoints)`
  should be greater than or equal to the largest event or follow-up time used in
  the model.

# Keywords

- `niter::Integer=1000`: Maximum number of iterations for the optimization
  algorithm.
- `tol::Float64=1e-4`: Convergence tolerance for the iterative estimation
  procedure.
- `factor::Float64=0.75`: Step-size adjustment factor used in the optimization
  routine. It should be between 0 and 1.

# Returns

A vector of named tuples with fields:

- `z_idx`: Covariate index.
- `beta`: Estimated regression coefficient.
- `se`: Estimated standard error.

# Details

This function fits a proportional hazards model for interval-censored
single-event data. The model assumes that the hazard function depends
multiplicatively on covariates through a log-linear form. The single-event
problem is represented internally as a two-state multistate model with one
transition, `1 => 2`.

# Examples

```julia
using IntervalCensoredMultistate

left = Float64[0.0, 0.2, 0.5, 0.4]
right = Float64[0.3, 0.6, 0.8, 1.0]
event = Bool[true, true, true, false]
zzi = Float64[
    0.0 1.0;
    1.0 0.0;
    0.5 1.0;
    1.0 1.0
]
timepoints = collect(Float64, 0.1:0.1:1.0)

summary = fit_single_event(left, right, event, zzi, timepoints; niter=100)
```
"""
function fit_single_event(
  left::Vector{Float64},
  right::Vector{Float64},
  event::Vector{Bool},
  zzi::Matrix{Float64},
  timepoints::Vector{Float64};
  niter::Integer=1000,tol::Float64=1e-4,factor::Float64=0.75)

  nind = length(left)
  length(right) == nind || throw(ArgumentError("`right` must have the same length as `left`."))
  length(event) == nind || throw(ArgumentError("`event` must have the same length as `left`."))
  size(zzi, 1) == nind || throw(ArgumentError("`zzi` must have one row per individual."))

  TTij = fill(NaN, nind, 3)
  ssij = fill(false, nind, 3, 2)
  nobss = Vector{Int64}(undef, nind)

  for ii in 1:nind
    if event[ii] && left[ii] > 0.0
      TTij[ii, :] = [0.0, left[ii], right[ii]]
      ssij[ii, 1, 1] = true
      ssij[ii, 2, 1] = true
      ssij[ii, 3, 2] = true
      nobss[ii] = 3
    elseif event[ii] && left[ii] == 0.0
      TTij[ii, 1:2] = [0.0, right[ii]]
      ssij[ii, 1, 1] = true
      ssij[ii, 2, 2] = true
      nobss[ii] = 2
    else
      TTij[ii, 1:2] = [0.0, left[ii]]
      ssij[ii, 1, 1] = true
      ssij[ii, 2, 1] = true
      nobss[ii] = 2
    end
  end

  F_zidx = [Int64[] for _ in 1:2, _ in 1:2]
  F_zidx[1, 2] = collect(Int64, 1:size(zzi, 2))
  possible_transition = Bool[false true; false false]

  result = fit_multistate(
    TTij, ssij, nobss, zzi, F_zidx, possible_transition, timepoints;
    niter=niter, tol=tol, factor=factor)

  return [(z_idx=row.z_idx, beta=row.beta, se=row.se) for row in result.summary]
end

"""
    fit_competing_risks(left, right, cause, zzi, F_zidx, timepoints; niter=1000, tol=1e-4, factor=0.75)

Fit a proportional hazards model for interval-censored competing-risks data.

Estimates regression coefficients in a proportional hazards framework for
interval-censored competing-risks data. Each individual is observed within an
interval `(left, right]`, with a recorded failure cause if an event occurs.
Cause-specific hazards are modeled using transition-specific covariate effects.

# Arguments

- `left::Vector{Float64}`: Left endpoints of observation intervals. For the
  `i`-th individual, `left[i]` is the lower bound of the interval in which the
  event time lies. In cases of right censoring, `left[i]` is the last follow-up
  time.
- `right::Vector{Float64}`: Right endpoints of observation intervals. For the
  `i`-th individual, `right[i]` is the upper bound of the interval in which the
  event time lies.
- `cause::Vector{Int64}`: Event type. `cause[i] == 0` indicates right
  censoring, while positive integers `1, 2, ..., K` indicate the cause of
  failure.
- `zzi::Matrix{Float64}`: Covariate matrix. The `i`-th row corresponds to the
  `i`-th individual and each column to a predictor variable.
- `F_zidx_by_cause::Vector{Vector{Int64}}`: Covariate indices for each cause.
  `F_zidx_by_cause[k]` contains the indices of covariates used to model the
  cause-specific hazard for cause `k`.
- `timepoints::Vector{Float64}`: Timepoints for the discrete approximation.
  Values must be finite, positive, and strictly increasing. `maximum(timepoints)`
  should be greater than or equal to the largest event or follow-up time used in
  the model.

# Keywords

- `niter::Integer=1000`: Maximum number of iterations for the optimization
  algorithm.
- `tol::Float64=1e-4`: Convergence tolerance for the iterative estimation
  procedure.
- `factor::Float64=0.75`: Step-size adjustment factor used in the optimization
  routine. It should be between 0 and 1.

# Returns

A vector of named tuples with fields:

- `cause`: Failure cause.
- `z_idx`: Covariate index.
- `beta`: Estimated regression coefficient.
- `se`: Estimated standard error.

# Details

This function fits cause-specific proportional hazards models under interval
censoring. For each failure cause, a separate hazard function is modeled with a
log-linear dependence on selected covariates. The competing-risks problem is
represented internally as a multistate model with transitions from state 1 to
one absorbing state per failure cause.

# Examples

```julia
using IntervalCensoredMultistate

left = Float64[0.0, 0.2, 0.5, 0.4]
right = Float64[0.3, 0.6, 0.8, 1.0]
cause = Int64[1, 2, 1, 0]
zzi = Float64[
    0.0 1.0;
    1.0 0.0;
    0.5 1.0;
    1.0 1.0
]
F_zidx_by_cause = Vector{Int64}[[1, 2], [1]]
timepoints = collect(Float64, 0.1:0.1:1.0)

summary = fit_competing_risks(
    left, right, cause, zzi, F_zidx_by_cause, timepoints;
    niter=100
)
```
"""
function fit_competing_risks(
  left::Vector{Float64},
  right::Vector{Float64},
  cause::Vector{Int64},
  zzi::Matrix{Float64},
  F_zidx_by_cause::Vector{Vector{Int64}},
  timepoints::Vector{Float64};
  niter::Integer=1000,tol::Float64=1e-4,factor::Float64=0.75)

  nind = length(left)
  length(right) == nind || throw(ArgumentError("`right` must have the same length as `left`."))
  length(cause) == nind || throw(ArgumentError("`cause` must have the same length as `left`."))
  size(zzi, 1) == nind || throw(ArgumentError("`zzi` must have one row per individual."))
  minimum(cause) >= 0 || throw(ArgumentError("`cause` values must be nonnegative."))

  nstates = maximum(cause) + 1
  length(F_zidx_by_cause) == nstates - 1 ||
    throw(ArgumentError("`F_zidx_by_cause` must contain one vector per positive cause."))

  TTij = fill(NaN, nind, 3)
  ssij = fill(false, nind, 3, nstates)
  nobss = Vector{Int64}(undef, nind)

  for ii in 1:nind
    if cause[ii] > 0 && left[ii] > 0.0
      TTij[ii, :] = [0.0, left[ii], right[ii]]
      ssij[ii, 1, 1] = true
      ssij[ii, 2, 1] = true
      ssij[ii, 3, cause[ii] + 1] = true
      nobss[ii] = 3
    elseif cause[ii] > 0 && left[ii] == 0.0
      TTij[ii, 1:2] = [0.0, right[ii]]
      ssij[ii, 1, 1] = true
      ssij[ii, 2, cause[ii] + 1] = true
      nobss[ii] = 2
    else
      TTij[ii, 1:2] = [0.0, left[ii]]
      ssij[ii, 1, 1] = true
      ssij[ii, 2, 1] = true
      nobss[ii] = 2
    end
  end

  F_zidx = [Int64[] for _ in 1:nstates, _ in 1:nstates]
  for cause_id in 1:(nstates - 1)
    F_zidx[1, cause_id + 1] = F_zidx_by_cause[cause_id]
  end

  possible_transition = fill(false, nstates, nstates)
  possible_transition[1, 2:nstates] .= true

  result = fit_multistate(
    TTij, ssij, nobss, zzi, F_zidx, possible_transition, timepoints;
    niter=niter, tol=tol, factor=factor)

  return [
    (cause=row.to - 1, z_idx=row.z_idx, beta=row.beta, se=row.se)
    for row in result.summary
  ]
end

"""
    fit_multistate(TTij, ssij, nobss, zzi, F_zidx, possible_transition, timepoints; niter=1000, tol=1e-4, factor=0.75)

Fit a proportional hazards model for interval-censored multistate data.

Estimates regression coefficients in a proportional hazards framework for
interval-censored multistate data. Individuals start in state 1 at time 0 and
are followed over time. At each observation time, the exact state may be
uniquely observed or only partially observed.

# Arguments

- `TTij::Matrix{Float64}`: Observation times. `TTij[i, j]` is the observation
  time for the `j`-th observation of the `i`-th individual. Only
  `TTij[i, 1:nobss[i]]` is used.
- `ssij::Array{Bool,3}`: Possible state occupation indicators.
  `ssij[i, j, s] == true` if the `i`-th individual may occupy state `s` at the
  `j`-th observation time.
- `nobss::Vector{Int64}`: Number of observation times for each individual.
  `nobss[i]` is the number of observations for the `i`-th individual.
- `zzi::Matrix{Float64}`: Covariate matrix. The `i`-th row corresponds to the
  `i`-th individual and each column to a predictor variable.
- `F_zidx::Matrix{Vector{Int64}}`: Transition-specific covariate indices. The
  `(s1, s2)` element contains the indices of covariates used to model the
  transition from state `s1` to state `s2`. Use `Int64[]` for transitions with
  no covariates or impossible transitions.
- `possible_transition::Matrix{Bool}`: Allowable state transitions.
  `possible_transition[s1, s2] == true` if and only if the transition from
  state `s1` to state `s2` is allowed.
- `timepoints::Vector{Float64}`: Timepoints for the discrete approximation.
  Values must be finite, positive, and strictly increasing. The wrapper maps
  entries of `TTij` to their nearest timepoint before calling `reg_est_bone`.

# Keywords

- `niter::Integer=1000`: Maximum number of iterations for the optimization
  algorithm.
- `tol::Float64=1e-4`: Convergence tolerance for the iterative estimation
  procedure.
- `factor::Float64=0.75`: Step-size adjustment factor used in the optimization
  routine. It should be between 0 and 1.

# Returns

A named tuple containing:

- `F_betaz_est`: Regression coefficient estimates for each transition.
- `F_II_betaz`: Information matrices for the regression coefficient estimates.
- `F_hazard0`: Baseline transition rates.
- `F_Eg`: Inferred state occupation probabilities.
- `F_eps`: Change in regression coefficient estimates at convergence.
- `intervals`: The time intervals defined by `timepoints`.
- `TTij_rounded`: Observation times rounded to the nearest timepoint.
- `TTij_int`: Integer timepoint indices passed to `reg_est_bone`.
- `summary`: A vector of named tuples with fields `transition`, `from`, `to`,
  `z_idx`, `beta`, and `se`.

# Details

We consider `N` individuals who are observed at possibly unevenly spaced time
points. For the `i`-th individual, the observation times are
`TTij[i, 1], ..., TTij[i, nobss[i]]`.

There are `ns` possible states, where `ns = size(possible_transition, 1)`. The
logical matrix `possible_transition` encodes allowable transitions:

```julia
possible_transition[s1, s2] == true
```

if and only if a transition from state `s1` to state `s2` is allowed.

At each observation time `TTij[i, j]`, the individual's state may be uniquely
determined or partially observed. If the state is uniquely known to be `s`, then
`ssij[i, j, s] == true` and all other components of `ssij[i, j, :]` are
`false`. If the individual may occupy multiple states, for example `s1` or
`s2`, then `ssij[i, j, s1]` and `ssij[i, j, s2]` are `true`, representing the
set of possible states at that time.

There are `nz` predictor variables, where `nz = size(zzi, 2)`.
Transition-specific covariate effects are specified through `F_zidx`, where
`F_zidx[s1, s2]` contains the indices of covariates used to model the transition
from state `s1` to state `s2`.

# Examples

```julia
using IntervalCensoredMultistate

TTij = Float64[
    0.0 0.3 0.8;
    0.0 0.5 NaN
]
ssij = fill(false, 2, 3, 3)
ssij[1, 1, 1] = true
ssij[1, 2, 1] = true
ssij[1, 3, 2] = true
ssij[2, 1, 1] = true
ssij[2, 2, 1] = true
nobss = Int64[3, 2]
zzi = Float64[
    0.0 1.0;
    1.0 0.0
]
F_zidx = [Int64[] for _ in 1:3, _ in 1:3]
F_zidx[1, 2] = Int64[1, 2]
possible_transition = Bool[
    false true  false;
    false false false;
    false false false
]
timepoints = collect(Float64, 0.1:0.1:1.0)

fit = fit_multistate(
    TTij, ssij, nobss, zzi, F_zidx, possible_transition, timepoints;
    niter=10
)
fit.summary
```
"""
function fit_multistate(
  TTij::Matrix{Float64},
  ssij::Array{Bool,3},
  nobss::Vector{Int64},
  zzi::Matrix{Float64},
  F_zidx::Matrix{Vector{Int64}},
  possible_transition::Matrix{Bool},
  timepoints::Vector{Float64};
  niter::Integer=1000,tol::Float64=1e-4,factor::Float64=0.75)

  validate_timepoints(timepoints)
  size(TTij, 1) == size(ssij, 1) == length(nobss) == size(zzi, 1) ||
    throw(ArgumentError("`TTij`, `ssij`, `nobss`, and `zzi` must describe the same number of individuals."))
  size(TTij, 2) == size(ssij, 2) ||
    throw(ArgumentError("`TTij` and `ssij` must have the same observation dimension."))
  size(possible_transition, 1) == size(possible_transition, 2) == size(ssij, 3) ||
    throw(ArgumentError("`possible_transition` must be square with one row and column per state."))
  size(F_zidx) == size(possible_transition) ||
    throw(ArgumentError("`F_zidx` must have the same dimensions as `possible_transition`."))

  TTij_int = fill(Int64(-1), size(TTij))
  TTij_rounded = fill(NaN, size(TTij))

  for ii in 1:size(TTij, 1)
    for jj in 1:nobss[ii]
      if TTij[ii, jj] == 0.0
        TTij_int[ii, jj] = 0
        TTij_rounded[ii, jj] = 0.0
      else
        idx = nearest_timepoint_index(TTij[ii, jj], timepoints)
        TTij_int[ii, jj] = idx
        TTij_rounded[ii, jj] = timepoints[idx]
      end
    end
  end

  F_betaz_est, F_II_betaz, F_hazard0, F_Eg, F_eps = reg_est_bone(
    TTij_int,
    ssij,
    nobss,
    zzi,
    F_zidx,
    possible_transition,
    length(timepoints);
    niter=niter, tol=tol, factor=factor)

  return (
    F_betaz_est=F_betaz_est,
    F_II_betaz=F_II_betaz,
    F_hazard0=F_hazard0,
    F_Eg=F_Eg,
    F_eps=F_eps,
    intervals=[(j=j, start=j == 1 ? 0.0 : timepoints[j - 1], stop=timepoints[j]) for j in eachindex(timepoints)],
    TTij_rounded=TTij_rounded,
    TTij_int=TTij_int,
    summary=summarize_fit(F_betaz_est, F_II_betaz, F_zidx, possible_transition))
end

function validate_timepoints(timepoints::Vector{Float64})
  isempty(timepoints) && throw(ArgumentError("`timepoints` must not be empty."))
  all(isfinite, timepoints) || throw(ArgumentError("`timepoints` must be finite."))
  all(timepoints .> 0.0) || throw(ArgumentError("`timepoints` must be positive."))
  all(diff(timepoints) .> 0.0) || throw(ArgumentError("`timepoints` must be strictly increasing."))
end

function nearest_timepoint_index(x::Float64, timepoints::Vector{Float64})
  isfinite(x) || throw(ArgumentError("Observation times must be finite, except unused entries may be NaN."))

  idx = searchsortedfirst(timepoints, x)
  if idx == 1
    return 1
  elseif idx > length(timepoints)
    return length(timepoints)
  elseif abs(x - timepoints[idx - 1]) <= abs(x - timepoints[idx])
    return idx - 1
  else
    return idx
  end
end

function summarize_fit(F_betaz_est, F_II_betaz, F_zidx::Matrix{Vector{Int64}}, possible_transition::Matrix{Bool})
  rows = NamedTuple[]

  for transition in findall(possible_transition)
    from, to = Tuple(transition)
    beta = F_betaz_est[from, to]
    se = sqrt.(LinearAlgebra.diag(LinearAlgebra.inv(-F_II_betaz[from, to])))

    for kk in eachindex(beta)
      push!(rows, (
        transition="($from, $to)",
        from=from,
        to=to,
        z_idx=F_zidx[from, to][kk],
        beta=beta[kk],
        se=se[kk]))
    end
  end

  return rows
end
