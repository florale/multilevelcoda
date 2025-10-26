# library(bayestestR)
# 
# data(sbp)
# cilr <- complr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# m <- brmcoda(complr = cilr,
#              formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#                                 wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#                               chain = 4, cores = 8, iter = 4000, seed = 123, 
#              save_pars = save_pars(all = TRUE))
# m0 <- brmcoda(complr = cilr,
#               formula = Stress ~ 1 + (1 | ID),
#               chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# m1 <- brmcoda(complr = cilr,
#              formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 + (1 | ID),
#              chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# m2 <- brmcoda(complr = cilr,
#               formula = Stress ~ wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#               chain = 4, cores = 8, iter = 4000, seed = 123, save_pars = save_pars(all = TRUE))
# 
# comparison <- bayesfactor_models(m$model, m1$model, m2$model, denominator = m0$model)
# as.matrix(comparison)
# comparison
