# library(bayestestR)
# 
# data(sbp)
# cilr <- complr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# m <- brmcoda(complr = cilr,
#              formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#                               chain = 4, cores = 8, iter = 4000, seed = 123, 
#              save_pars = save_pars(all = TRUE))
# m0 <- brmcoda(complr = cilr,
#               formula = Stress ~ 1 + (1 | ID),
#               chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# m1 <- brmcoda(complr = cilr,
#              formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 + (1 | ID),
#              chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# m2 <- brmcoda(complr = cilr,
#               formula = Stress ~ wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#               chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# comparison <- bayesfactor_models(m$model, m1$model, m2$model, denominator = m0$model)
# as.matrix(comparison)
# comparison
