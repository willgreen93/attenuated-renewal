# Attenuated Renewal Models for Epidemic Forecasting in Overdispersed Networks

This repository contains the code supporting the paper:

**Green WD, Overton C, Burdon M, Ward T, Funk S.**  
*Robust Forecasting of Infectious Disease Epidemics in MSM Networks Using Attenuated Renewal Processes*

---

## Overview

Forecasting emerging infectious disease outbreaks is particularly challenging in populations where transmission occurs on highly overdispersed contact networks, such as men who have sex with men (MSM). In these settings, early epidemic growth is often driven by a small number of highly connected individuals, leading to rapid early acceleration followed by structural attenuation as high-degree nodes are depleted.

This repository implements and evaluates a suite of short-term epidemic forecasting models, with a focus on attenuated renewal processes that explicitly account for this early depletion. Using both large-scale stochastic simulations on synthetic MSM-like networks and retrospective forecasts of the 2022 mpox outbreak in multiple cities, the code demonstrates that sigmoid-attenuated renewal models provide more stable and reliable forecasts than standard renewal or purely statistical approaches.

---

## Key features

- Simulation of epidemics on synthetic MSM-like contact networks with:
  - highly overdispersed degree distributions  
  - assortative mixing by degree and sexual role  
  - overlapping one-time and recurrent partnership networks
- Implementation of six forecasting models, spanning:
  - log-random walk baselines  
  - classical renewal models  
  - multiple forms of attenuation in the effective reproduction number \(R(t)\)
- Bayesian inference using **Hamiltonian Monte Carlo (Stan)**
- Systematic forecast evaluation using the **Continuous Ranked Probability Score (CRPS)**
- Application to real mpox incidence data from New York City, San Francisco, Madrid, and Barcelona

---

## Forecasting models implemented

The repository includes implementations of the following forecasting approaches:

- Log-random walk on incidence  
- Renewal model with random walk in R(t)  
- Biased random walk in R(t)
- Exponential attenuation of R(t) in time  
- Exponential attenuation of R(t) by attack rate  
- Power-law attenuation of R(t) by attack rate  
- Reverse-sigmoid attenuation of R(t) in time

---

## Setup

Before running the scripts, create the required output directories:

```r
dirs <- c("mixing_matrices", "simulation/sim_2", "simulation/sim_homo",
          "simulation/sim_het", "fits/sim/fit_sim3", "fits/sim/fit_homo",
          "fits/sim/fit_het", "fits/cities", "fits/outputs", "figures")
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
```

All scripts should be run from the project root directory.

---

## Data

The repository includes mpox incidence data from Barcelona and Madrid (`data/mpox_BCN_MAD.xlsx`), as well as case data from NYC and San Francisco.

---

## Session info

This code was developed using:

- R version 4.3.x
- Key packages: `cmdstanr`, `rstan`, `dplyr`, `ggplot2`, `igraph`, `EpiEstim`, `scoringutils`
