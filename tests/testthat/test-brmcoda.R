## Check prediction using a simulated dataset
nGrps <- 100
k <- 3
set.seed(1234)
u <- rnorm(nGrps, mean = 0, sd = 1)
day1.pa <- 2 + u + rnorm(nGrps, mean = 0, sd = .1)
day2.pa <- 3 + u + rnorm(nGrps, mean = 0, sd = .1)
day3.pa <- 4 + u + rnorm(nGrps, mean = 0, sd = .1)

day1 <- data.table(
  TST = 479, Wake = 961,
  Day = 1, PA = day1.pa, ID = seq_len(nGrps))
day2 <- data.table(
  TST = 480, Wake = 960,
  Day = 2, PA = day2.pa, ID = seq_len(nGrps))
day3 <- data.table(
  TST = 481, Wake = 959,
  Day = 3, PA = day3.pa, ID = seq_len(nGrps))
daydata <- rbind(day1, day2, day3)

sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = daydata, sbp = sbp,
                parts = c("TST", "Wake"), idvar = "ID")
psub <- possub(c("TST", "Wake"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = PA ~ wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
daydata2 <- cbind(daydata, fitted(m$Model))

test_that("wilr from brmcoda gives expected predictions", {
  expect_true(all.equal(daydata2[, PA[Day == 2] - PA[Day == 1]], 
                        daydata2[, Estimate[Day == 2] - Estimate[Day == 1]], tolerance = 0.2))
  expect_true(all.equal(daydata2[, PA[Day == 2] - PA[Day == 3]], 
                        daydata2[, Estimate[Day == 2] - Estimate[Day == 3]], tolerance = 0.2))
})