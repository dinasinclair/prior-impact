# Test file for Bayesian Cities 1D file

#library(testthat)
source('citiesBayesian1DwK.R')
 

context('testing if I can write tests')

test_that('test1', {
  expect_equal(I,3)
})

context('testing generate_basic')
# Test to see that the dimensions of everything are okay
# remember  basic_dat_generated <- list(I=I,Y=Y$mean,sigmaSq=Y$var)

test_that('Dimensions', {
  basic <- generate_basic(I=3,SF=1,mu=2,tau=2)
  expect_equal(basic$I,3)
  expect_equal(length(basic$Y), 3)
  expect_equal(length(basic$sigmaSq),3)
})

test_that('TauSq Zero', {
  basic <- generate_basic(I=3,SF=1,mu=2,tauSq=0)
  expect_equal(mean(basic$Y), 2)
})

test_that('SF (SigmaSq Mult Factor) Zero', {
  basic <- generate_basic(I=3,SF=0,mu=2,tauSq=2)
  expect_equal(mean(basic$sigmaSq), 0)
})

# Test that if var is 0, that we'd get what we'd expect both in the case where you change your mind
# the case where you wouldn't change your mind

# Test cases where you feed in new draws? Does that violate testing by testing the middle? I don't think so
# I think it should be okay to have a step post feeding it new draws and a step before

# Test that bayesian stuff is giving similar answers to a non-NUTS framework? Need to find one...
# This involves testing the stan file separately, tbh

# Test that we're getting all of the combinations

# Try setting the seed a few different ways, see if get the same results?

