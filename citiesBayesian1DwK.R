library("rstan")

generate_basic <- function(I,SF,mu,tauSq){
  
  # Generate theta and sigmaSq
  theta <- rnorm(I, mean = mu, sd = sqrt(tauSq)) # with theta ~ N(mu,tauSq)
  sigmaSq <- SF*runif(I) # with sigmaSq ~ U(0,1)
  Y <- list(mean=theta, var=sigmaSq)
  
  # Save our generated input data together in a list
  basic_dat_generated <- list(I=I,Y=Y$mean,sigmaSq=Y$var)
  
  # Display what we've generated
  return(basic_dat_generated)
}

extract_fit<- function(basic_dat_generated){
  Y <- list(mean = basic_dat_generated$Y, var = basic_dat_generated$sigmaSq)
  fit <- stan(file = 'randomEffectsModel1D.stan', 
              data = basic_dat_generated, 
              iter = 1000, chains = 2)
  
  # Readjust knowledge of Y based on REM
  params <- extract(fit)
  for (i in 1:length(Y$mean)){
    Y$mean[i] <- mean(params$theta[,i])
  }
  return(Y)
}

# Calculate new data for all K new pilots
get_pilot_results <- function(K,Y,Q){
  # Before update, Y_P is the same as Y
  Y_P <- Y
  
  # Update for each new pilot k
  for (k in K){
    # Gather New Pilot Data
    new_sigmaSq <- Y$var[k] * (1/Q)
    new_mean <- rnorm( 1 , mean = Y$mean[k] , sd = sqrt(new_sigmaSq ))
    
    # Combine the old and new data 
    post_pilot <- update_Y(Y$mean[k],new_mean,Y$var[k],new_sigmaSq)
    Y_P$mean[k] <- post_pilot$mean
    Y_P$var[k] <- post_pilot$var

  }
  return(Y_P)
}

update_Y <- function(mu1, mu2, sigSq1, sigSq2){
  update_mean  <- (mu1*sigSq2 + mu2*sigSq1)/(sigSq1 + sigSq2)
  update_var <- (sigSq1*sigSq2)/(sigSq1 + sigSq2)
  return(list(mean = update_mean , var = update_var))
}

new_ranking <- function(fit_updated,Y_P) {
  params_updated <- extract(fit_updated)
  Y_updated <- Y_P
  for (i in 1:length(Y_updated$mean)){
    Y_updated$mean[i] <- mean(params_updated$theta[,i])
  }
  new_rank <- order(Y_updated$mean, decreasing=TRUE)
  return(new_rank)
}

change_mind <- function(K,Y,original_rank,num_final_cities,I,Q){
  # This function returns TRUE if the pilots run in K change the final actions
  # and returns FALSE otherwise, given original ranking and data Y.
  
  # Simulate pilots
  Y_P <- get_pilot_results(K,Y,Q)
  # Update thetas
  updated_dat_generated <- list(I, Y=Y_P$mean, sigmaSq=Y_P$var)
  fit_updated <- stan(file = 'randomEffectsModel1D.stan', 
                      data = updated_dat_generated, 
                      iter = 1000, chains = 2)
  # Use new thetas to get a new city ranking
  new_rank <- new_ranking(fit_updated,Y_P) 
  # Take the top F cities from each ranking as final city choice
  original_choice <- original_rank[1:num_final_cities]
  new_choice <- new_rank[1:num_final_cities]
  # Return the setequality of the two rankings (set bc order doesn't matter)
  return(!setequal(original_choice,new_choice))
}

overall<-function(I, SF, num_pilots,num_final_cities,num_hypothetical_draws,Q,seed,mu,tau){
  
  set.seed(seed)
  basic_dat_generated <- generate_basic(I,SF,mu,tau)
  print(basic_dat_generated)
  Y <- extract_fit(basic_dat_generated)
  
  # Calculate the city ranking of the original data
  original_rank <- order(Y$mean, decreasing=TRUE)
  original_rank
  
  # Generate all combinations of cities $K$, store them in vector 'combinations'
  combinations <-combn(seq(I),num_pilots)
  print(combinations)
  # Initiate number of minds changed (nmc) to zero
  nmc <- numeric(ncol(combinations))
  
  # Loop through all combinations K, updating nmc num_hypothetical_draws times
  for (i in 1:ncol(combinations)){
    for (j in 1:num_hypothetical_draws){
      print("i is")
      print(i)
      K <- combinations[,i]
      nmc[i] <- nmc[i] + change_mind(K,Y,original_rank,num_final_cities)
    }
  } 
  
  print(Y)
  print("nmc")
  print(nmc)
  print(combinations)
  for (i in 1:length(nmc)){
    print(combinations[,i])
    print(paste("Number of times minds changed: ",nmc[i],"/",num_hypothetical_draws))
  } 
}

# data <- list(I=2,Y=c(2,2),sigmaSq=c(1,1))
# extract_fit(data)

# overall(I=3,
#         SF=1,
#         num_pilots=3,
#         num_final_cities=2,
#         num_hypothetical_draws=1,
#         Q=5,
#         seed=17,
#         mu=2,
#         tau=1)

# Y <- list(mean = c(0,0), var = c(1,1))
# Y_P <- get_pilot_results(K=c(1,2),Y=Y,Q=1)
# Y
# Y_P

Y <- list(mean = c(-1000,0,0,0), var = c(1,1,1,1))
change_mind(K=c(1),Y,original_rank=c(1,2,3,4),num_final_cities=1,I=4,Q=1)