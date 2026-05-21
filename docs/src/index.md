# IntervalCensoredMultistate.jl

IntervalCensoredMultistate.jl is a Julia package for fitting proportional
hazards models with interval-censored multistate data. It is designed for
settings where individuals are observed at discrete and possibly irregular time
points, and their true state occupations between observation times are not
directly observed.

The package provides high-level wrappers for common event-history models and a
lower-level computational core for custom multistate workflows:

- `fit_single_event` fits interval-censored single-event models.
- `fit_competing_risks` fits interval-censored competing-risks models.
- `fit_multistate` fits general interval-censored multistate models.
- `reg_est_bone` is the lower-level estimation routine used by the wrappers.

## Installation

```julia
using Pkg
Pkg.add("IntervalCensoredMultistate")
```

For development from a local checkout, use:

```julia
using Pkg
Pkg.develop(path="/path/to/IntervalCensoredMultistate")
```

## State Transition Data

We consider `N` individuals who all start in state 1 at time 0 and are followed
over time.

Each individual is observed at a sequence of possibly unevenly spaced time
points. For the `i`-th individual, the observation times are:

```julia
TTij[i, 1], ..., TTij[i, nobss[i]]
```

There are `ns` possible states in total. The allowable transitions between
states are specified by a user-defined transition matrix:

```julia
possible_transition[s1, s2] == true
```

if and only if a transition from state `s1` to state `s2` is allowed.

At each observation time `TTij[i, j]`, the individual's state may be either
uniquely determined or only partially observed. If the state is uniquely known
to be `s`, then `ssij[i, j, s] == true` and all other components of
`ssij[i, j, :]` are false. If the individual may occupy multiple states, for
example `s1` or `s2`, then `ssij[i, j, s1]` and `ssij[i, j, s2]` are true,
indicating a set of possible states.

## Covariates and Transition-Specific Effects

There are `nz` predictor variables from the `N` individuals used to model the
risks of state transitions. These covariates are stored in a matrix `zzi`, where
`zzi[i, k]` is the `k`-th covariate for the `i`-th individual.

Different sets of covariates can be used for different transitions. This is
specified through a matrix of vectors:

```julia
F_zidx[s1, s2]
```

The `(s1, s2)` entry contains the indices of covariates used to model the
transition from state `s1` to state `s2`. Use `Int64[]` for impossible
transitions or transitions with no covariates.

## Single-Event Example

For standard interval-censored time-to-event data, use `fit_single_event`.
Each event time is known to lie in `(left[i], right[i]]`; observations with
`event[i] == false` are treated as right-censored at `left[i]`.

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

The returned `summary` is a vector of named tuples with the covariate index
(`z_idx`), coefficient estimate (`beta`), and standard error (`se`).

## Competing-Risks Example

Competing risks are represented as transitions from the initial state to one
absorbing state per cause. A cause value of `0` indicates right censoring.

```julia
cause = Int64[1, 2, 1, 0]
F_zidx_by_cause = Vector{Int64}[[1, 2], [1]]

summary = fit_competing_risks(
    left, right, cause, zzi, F_zidx_by_cause, timepoints;
    niter=100
)
```

## Multistate Example

The general multistate wrapper accepts observation-time, state-possibility,
covariate, and transition-structure arrays directly.

```julia
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

## Reference

You, L., Liu, X., & Krischer, J. (2024). A discrete approximation method for
modeling interval-censored multistate data. *Statistics in Medicine*, 43(12),
2452-2471. PMID: 38599784, PMCID: PMC11109708, DOI: 10.1002/sim.10079.
