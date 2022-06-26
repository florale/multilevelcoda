# library(bayestestR)
# 
# data(sbp)
# cilr <- compilr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# m <- brmcoda(compilr = cilr,
#              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#                               chain = 4, cores = 8, iter = 4000, seed = 123, 
#              save_pars = save_pars(all = TRUE))
# m0 <- brmcoda(compilr = cilr,
#               formula = STRESS ~ 1 + (1 | ID),
#               chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# m1 <- brmcoda(compilr = cilr,
#              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + (1 | ID),
#              chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# m2 <- brmcoda(compilr = cilr,
#               formula = STRESS ~ wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#               chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# comparison <- bayesfactor_models(m$Model, m1$Model, m2$Model, denominator = m0$Model)
# as.matrix(comparison)
# comparison
