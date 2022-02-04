test_that("examplefun works", {
  x <- examplefun(1:100)

  ## check types
  expect_s3_class(x, "data.table")
  expect_type(x$x, "double")
})

test_that("examplefun errors if it should", {
  expect_error(example <- fun("does this fail?"))
})
