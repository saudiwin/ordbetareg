# run all scripts necessary to replicate & compile paper and SI

# whether to run the  simulations from scratch

run_sim <- F

if(run_sim)  {
  
  print("Now on random simulation")
  
  source("ordered_beta_reg_sim.R")
  
  print("Now on fixed simulation")
  
  source("ordered_beta_reg_sim_fixed.R")
  
}

print("Compiling paper")

rmarkdown::render("kubinec_ord_betareg_accepted.Rmd")

print("Compiling supplementary information")

rmarkdown::render("kubinec_ord_betareg_appendix_anon.Rmd")