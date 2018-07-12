# Test file for Bayesian Cities 1D file
source('citiesBayesian1DwK.R')
 

context('testing if I can write tests')

test_that('test1', {
  expect_equal(I,3)
})

context('testing generate_basic')

# Test to see that the dimensions of everything are okay
test_that('Dimensions', {
  basic <- generate_basic(I=3,SF=1,mu=2,tau=2)
  expect_equal(basic$I,3)
  expect_length(basic$Y, 3)
  expect_length(basic$sigmaSq,3)
})

# Test that if we set Tau Squared to 0, we get no variation in theta~(mu,0)
test_that('TauSq Zero', {
  basic <- generate_basic(I=3,SF=1,mu=2,tauSq=0)
  expect_equal(basic$Y, c(2,2,2))
})

# Test that if we say that sigma should be zero, that it is indeed zero
test_that('SF (SigmaSq Mult Factor) Zero', {
  basic <- generate_basic(I=3,SF=0,mu=2,tauSq=2)
  expect_equal(basic$sigmaSq, c(0,0,0))
})

# Test that bayesian stuff is giving similar answers to a non-NUTS framework? Try LAGS?
# This involves testing the stan file separately, tbh
context('testing stan fit and extraction')

test_that('Stan code generation', {
  # Fun fact: using 0 for sigma sq means it doesn't converge/work the way one would want. TODO look at why that would happen?
  data <- list(I=2,Y=c(2,2),sigmaSq=c(1,1))
  Y <- extract_fit(data)
  expect_true(mean(Y$mean) >= 1.5)
  expect_true(mean(Y$mean) <= 2.5)
  expect_equal(Y$var, c(1,1))
})

# Test that if var is 0, that we'd get what we'd expect both in the case where you change your mind
# the case where you wouldn't change your mind
context('testing change mind')

test_that('Expect to change mind if original idea is terrible', {
  Y <- list(mean = c(-1000,0,0,0), var = c(1,1,1,1))
  expect_true(change_mind(K=c(1),Y,original_rank=c(1,2,3,4),num_final_cities=1,I=4,Q=1))
})

test_that('Expect not to change mind if variance is really small # TODO this is a bad test, fix', {
  Y <- list(mean = c(400,300,200,100), var = c(.01,.01,.01,.01))
  expect_false(change_mind(K=c(1,2,3,4),Y,original_rank=c(1,2,3,4),num_final_cities=1,I=4,Q=1))
})

context('testing pilot generation')
# Test cases where you feed in new draws, check that we're updating the variance correctly.
test_that('Test K selection', {
  Y <- list(mean = c(0,0,0,0), var = c(1,1,1,1))
  Y_P <- get_pilot_results(K=c(1,3),Y,Q=1)
  expect_equal(Y_P$var, c(0.5,1,0.5,1))
  expect_false(isTRUE(all.equal(Y_P$mean, Y$mean)))
})

test_that('Test Q and variance function', {
  Y <- list(mean = c(0), var = c(1))
  Y_P <- get_pilot_results(K=c(1),Y,Q=10)
  expect_equal(Y_P$var, c(1/11))
  expect_false(isTRUE(all.equal(Y_P$mean, Y$mean)))
})

# Test that we're getting all of the combinations

# Try setting the seed a few different ways, see if get the same results?

