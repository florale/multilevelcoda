# test_that("bsub errors if it should", {
#   data(mcompd)
#   data(sbp)
# 
#   ps <- possub(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#   cilr <- compilr(data = mcompd, sbp = sbp, 
#                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#   
#   m <- brmcoda(compilr = cilr,
#                formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                core = 8, chain = 4)
#   
#   ## check types
#   x <- bsub(object = m, substitute = ps, minute = 5)
#   expect_type(x, "list")
#   expect_equal(length(x), length(m$CompIlr$parts))
# 
#   ## check errors for missing object
#   expect_error(x <- bsub(substitute = ps, minute = 5))
#   
#   ## check errors for missing substitute
#   expect_error(x <- bsub(objetc = m, minute = 5))
#   
#   ## check errors when reference grid has matching names with ILRs
#   rg <- data.table(bilr1 = 1)
#   expect_error(x <- bsub(object = m, substitute = ps, minute = 5, regrid = rg))
# })