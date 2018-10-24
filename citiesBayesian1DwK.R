library("rstan")
library("brms")

generate_basic <- function(I,mu,tau,SF=1){
  
  # Generate theta and sigmaSq
  theta <- rnorm(I, mean = mu, sd = tau) # with theta ~ N(mu,tau)
  sigmaSq <- SF*runif(I) # with sigmaSq ~ U(0,1)
  Y <- list(mean=theta, var=sigmaSq)
  
  # Save our generated input data together in a list
  basic_dat_generated <- list(I=I,Y=Y$mean,sigmaSq=Y$var)
  
  # Display what we've generated
  return(basic_dat_generated)
}

extract_fit<- function(data){
  # data has the form list(I=I,Y=Y$mean,sigmaSq=Y$var)
  Y <- list(mean = data$Y, var = data$sigmaSq)
  fit <- stan(file = 'randomEffectsModel1D.stan', 
              data = data, 
              iter = 1000, chains = 2, control=list(adapt_delta=0.80, max_treedepth=10))
  #pairs(fit)
  
  # Readjust knowledge of Y based on REM
  params <- extract(fit)
  for (i in 1:length(Y$mean)){
    Y$mean[i] <- mean(params$theta[,i])
  }
  return(list(Y=Y,mu=mean(params$mu), tau=mean(params$tau)))
}

extract_fit_brms<- function(data){
  # data has the form list(I=I,Y=Y$mean,sigmaSq=Y$var)
  Y <- list(mean = data$Y, var = data$sigmaSq)
  
  # Okay now the data does not have to be in the weird list form. we do need
  # Y_city ~ (mu|city) + (error|city) + population-wide effects
  # mu for each, group(city) for each, error for each, Y city as output
  fit <- brm(Y_new ~ Y|mu, 
              data = data, 
              family = lognormal)
  #pairs(fit)
  
  # Readjust knowledge of Y based on REM
  params <- extract(fit)
  for (i in 1:length(Y$mean)){
    Y$mean[i] <- mean(params$theta[,i])
  }
  return(list(Y=Y,mu=mean(params$mu), tau=mean(params$tau)))
}

# Calculate new data for all K new pilots
get_pilot_results <- function(K,Y,Q=1){
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

change_mind <- function(K,Y,original_rank,num_final_cities,Q=1){
  # This function returns TRUE if the pilots run in K change the final actions
  # and returns FALSE otherwise, given original ranking and data Y.
  
  # Simulate pilots
  Y_P <- get_pilot_results(K,Y,Q)
  
  # Update thetas
  updated_data <- list(I=length(Y$mean), Y=Y_P$mean, sigmaSq=Y_P$var)
  Y_updated <- extract_fit(updated_data)$Y
  
  # Use new thetas to get a new city ranking
  new_rank <- order(Y_updated$mean, decreasing=TRUE)
  
  # Take the top F cities from each ranking as final city choice
  original_choice <- original_rank[1:num_final_cities]
  new_choice <- new_rank[1:num_final_cities]
  
  # Return the setequality of the two rankings (set bc order doesn't matter)
  return(!setequal(original_choice,new_choice))
}

overall<-function(data,num_pilots,num_final_cities,num_draws,Q=1){

  # Do initial REM fit
  Y <- extract_fit(data)$Y
  
  # Calculate the city ranking of the original data
  original_rank <- order(Y$mean, decreasing=TRUE)
  
  # Generate all combinations of cities $K$, store them in vector 'combinations'
  combinations <-combn(seq(data$I),num_pilots)

  # Initiate number of minds changed (nmc) to zero
  nmc <- numeric(ncol(combinations))
  
  # Loop through all combinations K, updating nmc num_draws times
  for (i in 1:ncol(combinations)){
    for (j in 1:num_draws){
      K <- combinations[,i]
      nmc[i] <- nmc[i] + change_mind(K,Y,original_rank,num_final_cities)
    }
  } 
  
  return(list(nmc=nmc, combinations=combinations))
}

#generate_basic(I=3,mu=0,tau=1,SF=1)

