test_that("compilr errors if it should", {
  expect_error(example <- compilr())
})