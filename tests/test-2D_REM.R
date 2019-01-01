# Test file for Bayesian Cities 1D file
library("metafor")
source("C:/Users/Dina/Documents/Oxford/Evaluating Prior Impact/evaluating_2D_REM.Rmd")
set.seed(17)

context('testing generate_basic')

# Test to see that the dimensions of everything are okay
test_that('Dimensions', {
  basic <- generate_grouped_basic(I=3,mu=2,tau=2,N=4)
  expect_equal(basic$I,3)
  expect_equal(basic$N,4)
  expect_length(basic$Y_mean, 3)
  expect_length(basic$Y_sd,3)
  expect_length(basic$groups,3)
})

# Test that if we set all variances to 0 that we get the same mean in all Y
test_that('All Variance 0', {
  basic <- generate_basic(I=3,N=1,mu=2,tau=0, SF=0)
  expect_equal(basic$Y_mean, c(2,2,2))
})

# Test that if we set only tau to 0 but not group variance that we don't get the same mean in all Y
test_that('Tau 0 but group level still has variance', {
  basic <- generate_basic(I=3, N=1, mu=2,tau=0, SF=1)
  expect_not_equal(basic$Y_mean, c(2,2,2))
})

context('testing stan fit and extraction')

test_that('Only one data point', {
  data <- list(I=1, N=1, Y_mean=0, Y_sd=1, groups=1)
  fit <- extract_fit(data)
  expect_true(mean(fit$Y_mean) >= -0.5)
  expect_true(mean(fit$Y_mean) <= 0.5)
})

test_that('One group, two data points', {
  data <- list(I=2, N=1, Y_mean=c(0,0), Y_sd=c(1,1), groups=c(1,1))
  fit <- extract_fit(data)
  expect_true(mean(fit$Y_mean) >= -0.5)
  expect_true(mean(fit$Y_mean) <= 0.5)
})