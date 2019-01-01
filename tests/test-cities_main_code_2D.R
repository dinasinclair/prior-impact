# Test file for Bayesian Cities 2D file
library("metafor")
source("C:/Users/Dina/Documents/Oxford/Evaluating Prior Impact/cities_main_code_2D.R")
set.seed(17)

context('testing generate_basic')

# Test to see that the dimensions of everything are okay
test_that('Dimensions', {
  basic <- generate_basic(I=3,SF=1,mu=2,tau=2)
  expect_equal(basic$I,3)
  expect_length(basic$Y, 3)
  expect_length(basic$sigmaSq,3)
})

# Test that if we set Tau Squared to 0, we get no variation in theta~(mu,0)
test_that('Tau Zero', {
  basic <- generate_basic(I=3,SF=1,mu=2,tau=0)
  expect_equal(basic$Y, c(2,2,2))
})

context('testing extract_fit')

# test if we get the same result as before if we set G = 1
test_that('One Group is like a one dimensional model', {
  testData <- list(
    I = 3,
    G = 1,
    Y = c(1,2,3),
    sigmaSq = c(1,1,1),
    groups = c(1,1,1)
  )
})

# test if we have the exact same data for all points within one group
# that we get that result for that group
test_that('All group points identical leads to expected result', {
  testData <- list(
    I = 4,
    G = 2,
    Y = c(1,1,2,2),
    sigmaSq = c(0,0,0,0),
    groups = c(1,1,2,2)
  )
})