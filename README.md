README
================
Robert Kubinec
March 1st, 2020

This repository contains data and code for the paper, “Ordered Beta
Regression: A Parsimonious, Well-Fitting Model for Survey Sliders and
Visual Analog Scales.” The repository includes the following files:

  - `kubinec_ord_betareg.Rmd` A reproducible Rmarkdown document that can
    be run in Rstudio to re-produce the results. Note that the `data`
    folder in this repository contains necessary data to reproduce
    results.
  - `estimate_with_brms.Rmd` and `estimate_with_brms.html` These files
    show how to run ordered beta regression models using the R package
    [`brms`](https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf).
  - `define_ord_betareg.R` This R script contains all the auxiliary code
    needed to fit the model with R package `brms` (see vignette above
    for more info).
  - `ordered_beta_reg_sim.R` This R script will run a simulation
    comparing the ordered beta regression model to alternatives,
    including the zero-one-inflated Beta regression model (ZOIB). The
    output of a 10,000 run of this simulation is in `data/`.
  - `beta_logit.stan` This file contains the Stan code used to fit an
    ordered beta regression model in Stan.
  - `zoib.stan` This file contains Stan code used to fit the
    zero-one-inflated beta regression model (ZOIB).
  - `beta_logit_phireg.stan` This file constains Stan code to fit an
    ordered beta regression model with additional predictors for phi,
    the scale parameter in the distribution. These additional parameters
    allow for understanding the effect of covariates on encouraging
    clustered or more dispersed (estreme) responses from respondents.
  - `beta_logit_infl*.stan` These additional Stan files are various ways
    of parameterizing the midpoint of the scale when the midpoint is
    considered missing data. None of them appear to do a better job at
    predicting the outcome than versions that considered the midpoint to
    be observed data.
