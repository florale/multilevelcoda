## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "CairoPNG", dpi = 150, fig.path = "mlmcoda-"
)

## ----setup--------------------------------------------------------------------
library(multilevelcoda)
library(brms)
library(bayestestR)

options(digits = 3)

## ----data---------------------------------------------------------------------
data("mcompd") 
data("sbp") 

## -----------------------------------------------------------------------------
cilr <- compilr(data = mcompd, sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

## -----------------------------------------------------------------------------
head(cilr$TotalILR)

## ---- results = "hide"--------------------------------------------------------
mv <- brmcoda(compilr = cilr,
              formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ STRESS + (1 | ID),
              cores = 8, seed = 123, backend = "cmdstanr")

## -----------------------------------------------------------------------------
summary(mv$Model)

## -----------------------------------------------------------------------------
# intercept only
mv0 <- brmcoda(compilr = cilr,
               formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ 1 + (1 | ID),
               iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
               backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# full model
mv <- brmcoda(compilr = cilr,
              formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ STRESS + (1 | ID),
              iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
              backend = "cmdstanr", save_pars = save_pars(all = TRUE))

## -----------------------------------------------------------------------------
bayesfactor_models(mv$Model, denominator = mv0$Model)

