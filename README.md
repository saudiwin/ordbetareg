README
================
Robert Kubinec
January  11th, 2022

This repository contains data and code for the paper, “Ordered Beta Regression: A Parsimonious, Well-Fitting Model for Continuous Data with Lower and Upper Bounds", which is  now forthcoming at the journal *Political Analysis*. An ungated preprint can be found here: https://osf.io/preprints/socarxiv/2sx6y/. Replication files can be found both on Dataverse and Github. 

To replicate the paper, please first run the `install.R` script to make sure all relevant packages are installed. The script will also install `cmdstanr` and a version of `cmdstan`, which is the underlying MCMC sampling library from the Stan project. Installing `cmdstan` requires the R toolchain; if you have any trouble or are unsure see the `cmdstanr` package installation instructions: https://mc-stan.org/cmdstanr/articles/cmdstanr.html.

The file `master.R` will then run all the necessary scripts to compile the paper and supplementary information  (compilation requires a working Latex installation). Note that `master.R` by default loads the existing simulation data in the `data` folder. To fully reproduce the simulation, set the `run_sim` variable in `master.R` to `TRUE`. Note that running the full simulation can require up to a few days on a machine with ~40 cores. 

The repository includes the following files:

  - `kubinec_ord_betareg_accepted.Rmd` The accepted version of the reproducible Rmarkdown document that can
    be run in Rstudio to re-produce the results. Note that the `data`
    folder in this repository contains necessary data to reproduce
    results.
  - `estimate_with_brms.Rmd` and `estimate_with_brms.html` These files
    show how to run ordered beta regression models using the R package
    [`brms`](https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf).
  - `define_ord_betareg.R` This R script contains all the auxiliary code
    needed to fit the model with R package `brms` (see vignette above
    for more info).
  - `*_fit.rds` Fitted model object files to reproduce paper results much faster
  - `data/sim_cont_X*.RData` Simulation results to reproduce paper results much faster
  - `data/suffrage_paper_replicationfiles/EER-D-13-00718R2_maindata_suffrage.dta` Data from Toke and Aidt (2012)
  - `ordered_beta_reg_sim.R` This R script will run a simulation
    comparing the ordered beta regression model to alternatives,
    including the zero-one-inflated Beta regression model (ZOIB). The
    output of a 10,000 run of this simulation is saved in `data/` as `sim_cont_X.RData`.
  - `ordered_beta_reg_sim_fixed.R` This R script will run a simulation
    comparing the ordered beta regression model to alternatives, but with fixed rather than random draws 
    of relevant parameters (results are in the SI, not main paper). The
    output of a 4,000 run of this simulation is saved in `data/` as `sim_cont_X_fixed.RData`.
  - `beta_logit.stan` This file contains the Stan code used to fit an
    ordered beta regression model in Stan.
  - `zoib_nophireg.stan` This file contains Stan code used to fit the
    zero-one-inflated beta regression model (ZOIB).
  - `beta_logit_phireg.stan` This file constains Stan code to fit an
    ordered beta regression model with additional predictors for phi,
    the scale parameter in the distribution. These additional parameters
    allow for understanding the effect of covariates on encouraging
    clustered or more dispersed (estreme) responses from respondents.
  - `frac_logit.stan` This file contains a Stan parameterization of the 
    fractional logit model.
  - `beta_logit_infl*.stan` These additional Stan files are various ways
    of parameterizing the midpoint of the scale when the midpoint is
    considered missing data. None of them appear to do a better job at
    predicting the outcome than versions that considered the midpoint to
    be observed data.
  - `BibTexDatabase.bib` References necessary  to compile the paper.
  - `preamble.tex` Latex packages for paper
