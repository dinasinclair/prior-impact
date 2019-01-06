# This file creates data for a two-layer fit with one grouping level.
# It generates fake data, extracts a fit from that data using a stan REM model.
# TODO(dsinc): separate into different files
# TODO(dsinc): add tests
# TODO(dsinc): decide on sd vs tau conventions plus make sure we're not mixing sd vs var


library("rstan")

generateGroupedData <- function(I, mu, sd, N, SF=1){
  # Generates example dataset where data is grouped into N groups.
  # Note: Due to random group assignment some groups may have no data points assigned to them.
  # 
  # Arguments:
  #   I: the number of data points (cities)
  #   mu: the hyperparameter mean
  #   sd: the hyperparameter standard deviation
  #   N: the number of groups
  #   SF: the sigma factor, how wide spread we expect sigma to be

  # Generate mu/sd for each group
  mu_groups <- rnorm(N, mean = mu, sd = sd) # with theta ~ N(mu,sd)
  sd_groups <- SF * runif(N) # TODO(dsinc): correct sd generation to mirror stan code (lognormal?)
  G <- list(mu = mu_groups, sd=sd_groups)
  
  # Generate mu/sd for each study under a group
  group_assignment <- sample(1:N, I, replace=T)
  mu_studies <- numeric(I)
  sd_studies <- SF*runif(I)
  for (i in 1:I){
    mu_studies[i] <- rnorm(1,mean = G$mu[group_assignment[i]], sd = G$sd[group_assignment[i]])
  }
  Y <- list(mu = mu_studies, sd = sd_studies)
  
  # Save our generated input data together in a list
  generated_data <- list(I = I, N = N, Y = Y, G = G, groups = group_assignment)
  
  # Return what we've generated
  return(generated_data)
}

parseGeneratedGroupedData <- function(generated_data){
  # Given a data input, outputs the needed format for the stan file
  # Save our generated input data together in a list
  stan_input_data_format <- list(I = generated_data$I,
                                 N = generated_data$N,
                                 Y_mean = generated_data$Y$mu,
                                 Y_sd = generated_data$Y$sd,
                                 groups=generated_data$groups)
  return(stan_input_data_format)
}

extract_fit <- function(data){
  # This function returns the parameters backed out from the stan REM model. 
  # TODO(dsinc): should this instead be a function that takes in a bunch of arguments rather than one list?
  #
  # Arguments:
  #   data: a list in the form list(I, N, Y_mean, Y_sd, groups)
  # 
  # Returns:
  #   a list of updated parameters Y_mean Y_sd G_mean G_sd mu tau
  
  fit <- stan(file = 'randomEffectsModel2D.stan', 
              data = data, 
              iter = 1000, chains = 2, control=list(adapt_delta=0.99, max_treedepth=10))
  
  # Readjust knowledge of Y based on REM
  params <- extract(fit)
  Y <- list(mean = data$Y_mean, sd = data$Y_sd)
  G<- list(mean = data$G_mean, sd = data$G_sd)
  
  for (i in 1:length(data$Y_mean)){
    Y$mean[i] <- mean(params$theta[,i])
  }
  for (n in 1:length(data$groups)){
    G$mean[n] <- mean(params$G_mean[,n])
    G$sd[n] <- mean(params$G_sd[,n])
  }
  
  return(list(Y = Y, G = G, mu = mean(params$mu), tau = mean(params$tau)))
}


get_pilot_results <- function(K,Y,Q=1){
  # Calculate new data for all K new pilots
  # 
  # Arguments:
  #   K: int number of pilot studies to perform
  #   Y: (mu,sd) of all studies post first bayesian update
  #   Q: Factor by which pilot data uncertainty is lower than original data uncertainty
    
  # Before update, initialize Y_P (the post-pilot results) to Y
  Y_P <- Y
  
  # Update for each new pilot k
  for (k in K){
    # Gather new pilot data
    new_sd <- Y$sd[k] * (1/Q)
    new_mean <- rnorm( 1, mean = Y$mean[k], sd = new_sd)
    
    # Combine the old and new data 
    post_pilot <- update_Y(Y$mean[k], new_mean, Y$sd[k], new_sd)
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

data <- parseGeneratedGroupedData(generateGroupedData(I=10,mu=0,sd=10,N=3, SF=1))
fit <- stan(file = 'randomEffectsModel2D.stan', 
            data = data, 
            iter = 1000, chains = 2)

extract_fit(data)
