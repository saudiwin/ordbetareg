# This file will install all packages and also run cmdstanr's install_cmdstan function if it cannot 
# detect an existing cmdstan installation

# code borrowed from Vikram Baliga https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/

## First specify the packages of interest
packages = c("dplyr","rstanarm","tidyr",
             "lubridate","loo","kableExtra",
             "bayesplot","patchwork","stringr","grDevices","emojifont",
             "latex2exp","haven","ggplot2","moments",
             "posterior","brms","remotes","future.apply",
             "faux","rmarkdown","bookdown","tinytex","extrafont","binom","Hmisc",
             "ggthemes","ggtext")

print("Checking and installing packages.")

# create directory for data

dir.create("data")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    print(paste0("Now checking ", x))
    if (!require(x, character.only = TRUE)) {
      print(paste0("Package ",x," not installed, installing from CRAN."))
      install.packages(x, dependencies = TRUE)
    }
  }
)

# install cmdstanr

if (!require("cmdstanr", character.only = TRUE)) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}

# check for cmdstan installation

library(cmdstanr)

# checking if path to cmdstan exists
check_path <- try(cmdstan_path())

if('try-error' %in% class(check_path)) {
  
  user_input <- menu(c("Yes", "No"), title="Cmdstanr cannot find a cmdstan installation (if it does exist, you can also pass it as an option to set_cmdstan_path()). Would you like to install cmdstan from source? This is required to replicate the paper results.")
  
}

if(user_input==1L) {
  
  print("Installing cmsdstan.")
  
  print("If this fails, see package documentation at https://mc-stan.org/cmdstanr/.")
  
  install_cmdstan()
}



