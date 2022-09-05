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
library(doFuture)

options(digits = 3) # reduce number of digits shown

## ----data---------------------------------------------------------------------
data("mcompd") 
data("sbp")
data("psub")

## -----------------------------------------------------------------------------
cilr <- compilr(data = mcompd, sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

## ---- results = "hide"--------------------------------------------------------
m <- brmcoda(compilr = cilr,
             formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
               wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
             cores = 8, seed = 123, backend = "cmdstanr")

## -----------------------------------------------------------------------------
summary(m$Model)

## ---- results = "hide"--------------------------------------------------------
# intercept only model
m0 <- brmcoda(compilr = cilr,
             formula = STRESS ~ 1 + (1 | ID),
             iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
             backend = "cmdstanr", save_pars = save_pars(all = TRUE))

# between-person composition only model
m1 <- brmcoda(compilr = cilr,
             formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + (1 | ID),
             iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
             backend = "cmdstanr", save_pars = save_pars(all = TRUE))

# within-person composition only model
m2 <- brmcoda(compilr = cilr,
             formula = STRESS ~ wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
             iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
             backend = "cmdstanr", save_pars = save_pars(all = TRUE))

# full model
m <- brmcoda(compilr = cilr,
             formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
               wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
             iter = 6000, chains = 8, cores = 8, seed = 123, warmup = 1000,
             backend = "cmdstanr", save_pars = save_pars(all = TRUE))

## -----------------------------------------------------------------------------
comparison <- bayesfactor_models(m$Model, m1$Model, m2$Model, denominator = m0$Model)
comparison

## -----------------------------------------------------------------------------
update(comparison, reference = 1)

## -----------------------------------------------------------------------------
bsubm <- bsub(objec = m, substitute = psub, minute = 5)

## ---- results = "asis"--------------------------------------------------------
knitr::kable(bsubm$TST[abs(MinSubstituted) == 5])

## ----mlmcoda-plotbsubm2, fig.width = 9, fig.height = 6------------------------
plotsub(data = bsubm$TST, x = "sleep", y = "stress")

## -----------------------------------------------------------------------------
# Within-person substitution
wsubm <- wsub(objec = m, substitute = psub, minute = 5, summary = TRUE)

## ---- results = "asis"--------------------------------------------------------
knitr::kable(wsubm$TST[abs(MinSubstituted) == 5])

## ----mlmcoda-plotwsubm2, fig.width = 9, fig.height = 6------------------------
plotsub(data = wsubm$TST, x = "sleep", y = "stress")

