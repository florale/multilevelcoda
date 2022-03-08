test_that("compilr errors if it should", {
  expect_error(example <- compilr())
})

## test_that("compilr works", {
##   data(mcompd)
##   ## TODO need to create / define sbp
##   x <- compilr(data = mcompd[, 1:6], sbp = sbp, idvar = "ID")
##   ## check types
##   expect_type(x, "list")
##   expect_equal(length(x), 6L)
## })
