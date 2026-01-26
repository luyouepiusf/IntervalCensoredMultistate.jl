# IntervalCensoredMultistate.jl

IntervalCensoredMultistate is a Julia package for fitting proportional hazards models with interval-censored multistate data. The package is designed for settings where individuals are observed at discrete and possibly irregular time points, and their true state occupations between observation times are not directly observed.

## Key Features

* **Flexible State Space:** Define any transition structure (e.g., illness-death, recovery, or competing risks).
* **Transition-Specific Covariates:** Map specific predictors to specific state transitions.
* **Interval Censoring Support:** Handles data where state occupancy is only known at discrete, potentially unevenly spaced observation times.
* **State Occupancy Uncertainty:** Supports cases where the exact state at a given observation time is unknown. Users can specify a set of multiple "possible" states for each observation, allowing for more flexible interval censoring.

## Introduction of the problem

### State transition data

We consider `N` individuals who all start in state 1 at time 0 and are followed up to time `T` \(T = `ntimepoints`\).

Each individual is observed at a sequence of possibly unevenly spaced time points. For the `i`-th individual, the observation times are:

`TTij[i, 1]`, ..., `TTij[i, nobss[i]]`

There are **ns possible states** in total. The allowable transitions between states are specified by a user-defined transition matrix:

`possible_transition[s1, s2] = true`

if and only if a transition from state s1 to state s2 is allowed.

At each observation time `TTij[i, j]`, the individual's state may be either uniquely determined or only partially observed. If the state is uniquely known to be `s`, then `ssij[i, j, s]` is true and all other components of `ssij[i, j, :]` are false. If the individual may occupy multiple states (e.g., `s1` or `s2`), then `ssij[i, j, s1]` and `ssij[i, j, s2]` are true, indicating a set of possible states.

### Covariates and transition-specific effects

There are `nz` predictor variables (covariates) from the `N` individuals used to model the risks of state transitions. These covariates are stored in a matrix `zzi`, where:

`zzi[i, k]` is the `k`-th covariate for the `i`-th individual.

Importantly, the package allows **different sets of covariates for different transitions**. This is specified through a matrix of index vectors:

`F_zidx[s1, s2]` ⊆ {1, ..., `nz`}

which indicates which covariates are used to model the transition from state `s1` to state `s2`.

This structure allows highly flexible transition-specific regression models.

## Design philosophy

This package is intentionally minimal. The package currently exposes only one main function. It does not aim to provide a full modeling ecosystem with plotting, diagnostics, or simulation utilities. Instead, it provides a single, well-defined computational core that can be called from custom workflows, simulation studies, or higher-level modeling frameworks. Users are expected to handle data preprocessing, result interpretation, and downstream analysis externally.

## Applications

This framework is suitable for applications such as disease progression modeling, longitudinal clinical studies, reliability analysis, and other settings involving partially observed state processes.

## Installation

```julia
using Pkg
Pkg.add("IntervalCensoredMultistate") # Replace with actual name or local path
```

## Usage Guide

1. Define the Transition Logic
Use a BitMatrix to define which transitions are allowed between states.

```julia
# Example: 3-state model (1: Healthy, 2: Sick, 3: Dead)
# Transitions: 1->2, 1->3, and 2->3
possible_transition = [
    false  true   true;
    false  false  true;
    false  false  false
]
```

2. Map Predictors to Transitions
F_zidx allows you to specify which columns of your covariate matrix (zzi) influence which transitions.

```julia
# Assume zzi has 2 columns: [Age, Treatment]
F_zidx = Matrix{Vector{Int}}(undef, 3, 3)

F_zidx[1, 2] = [1, 2] # 1->2 influenced by Age and Treatment
F_zidx[1, 3] = [1]    # 1->3 influenced only by Age
F_zidx[2, 3] = [1]    # 2->3 influenced only by Age
```

3. Estimate Coefficients
Call the core estimation function with your prepared matrices.

results = estimate_coefficients(
    TTij,                # Observation times matrix (N x MaxObs)
    ssij,                # State possibility indicator (N x MaxObs x ns)
    nobss,               # Vector of observation counts per individual
    zzi,                 # Covariate matrix (N x nz)
    F_zidx,              # Mapping of zzi to transitions
    possible_transition, # The transition logic matrix
    ntimepoints          # Total time horizon (T)
)

---