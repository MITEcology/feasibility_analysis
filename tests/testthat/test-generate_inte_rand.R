test_that("output is a S x S matrix", {
  expect_equal(dim(generate_inte_rand(4, 1, 1, "norm")), 
               c(4,4))
})

test_that("output is a matrix", {
  expect_true(is.matrix(generate_inte_rand(4, 1, 1, "norm")))
})

test_that("output is a numeric", {
  expect_true(is.numeric(generate_inte_rand(4, 1, 1, "norm")))
})
