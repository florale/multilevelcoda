
# Packages
library(testthat)
library(data.table)
library(multilevelcoda)
library(extraoperators)
library(brms)
library(lme4)

# Model
#---------------------------------------------------------------------------------------------------
data(mcompd)
data(sbp)
data(psub)

cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")


test_that("complr errors if it should", {
  
  ## incorrect dataset type
  d <- list(
    bilr = 1:10,
    wilr = 1:10
  )
  expect_error(x <- complr(d, sbp = sbp, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID" ))
  
  ## incorrect sbp
  sbp1 <- matrix(c(
    1, -1, -1,-1,
    0, 1, -1, -1,
    0, 0, 1, -1), ncol = 4, byrow = TRUE)
  
  expect_error(x <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
                          sbp = sbp1, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID" ))
  
  ## incorrect sbp type
  sbp2 <- data.table(v1 = c(1, -1, -1,-1),
                     v2 = c(0, 1, -1, -1),
                     v3 = c(0, 0, 1, -1))
  
  expect_error(x <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
                          sbp = sbp2, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID" ))
  
  ## dataset contains ilr
  mcompd[, bilr1 := .1]
  expect_error(x <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
                           sbp = sbp, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID" ))
  
  # sbp contains invalid values
  sbpe <- matrix(c(
    1, 1, -1,-1, -1,
    1, -1, 0, 0, 0,
    0, 0, 1, -1, -2,
    0, 0, 0, 1, -1), ncol = 5, byrow = TRUE)
  expect_error(x <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
                           sbp = sbpe, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID"))
  
  # data has 0
  mcompd0 <- mcompd[ID %in% 1:10, .SD[1:3]][, .(TST = 0, WAKE = 0, MVPA = 0, LPA = 0, SB = 0)]
  expect_error(x <- complr(data =mcompd0,
                           sbp = sbp, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID"))
  })
