# Test file for Bayesian Cities 1D file

library(testthat)
source('citiesBayesian1DwK.Rmd')
 

context('testing if I can write tests')

test_that('test1', {
  expect_equal(I,3)
})


