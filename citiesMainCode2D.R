# This file creates data for a two-layer fit with one grouping level.
# It generates fake data, extracts a fit from that data using a stan REM model.

# DOCUMENTATION/CODE DEBT
# TODO(dsinc): separate into different files
# TODO(dsinc): add tests
# TODO(dsinc): decide on sd vs tau conventions plus make sure we're not mixing sd vs var
# TODO(dsinc): call it K or call it num_final_cities? (K can mean the city set as opposed to the number...)

# BUGS/CODE ERRORS
# TODO(dsinc): for some reason R is recompiling every time it runs the stan model, which makes it WAY too slow.


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

extractFit <- function(data){
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
  
  # Update values of Y, G, mu, tau appropriately
  for (i in 1:length(data$Y_mean)){
    Y$mean[i] <- mean(params$theta[,i])
  }
  for (n in 1:length(data$G_mean)){
    G$mean[n] <- mean(params$G_mean[,n])
    G$sd[n] <- mean(params$G_sd[,n])
  }
  mu = mean(params$mu)
  tau = mean(params$tau)
  
  return(list(Y = Y, G = G, mu = mu, tau = tau))
}


getPilotResults <- function(K,Y,Q=1){
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
    post_pilot <- combineMeanAndVar(mu1 = Y$mean[k], mu2 = new_mean, var1 = (Y$sd[k])^2, var2 = new_sd^2)
    Y_P$mean[k] <- post_pilot$mean
    Y_P$sd[k] <- sqrt(post_pilot$var)
    
  }
  return(Y_P)
}

combineMeanAndVar <- function(mu1, mu2, var1, var2){
  # Combines the means and variances of two data points (mu1, var1) (mu2, var2)
  
  combined_mean  <- (mu1*var2 + mu2*var1)/(var1 + var2)
  combined_var <- (var1*var2)/(var1 + var2)
  
  return(list(mean = combined_mean , var = combined_var))
}

# getNewRanking <- function(fit_updated,Y_P) {
#   # Returns a new ranking of programs by quality in decreasing order
#   # Arguments:
#   #   fit_updated: fit after being updated by TODO(dsinc): pilot? something? double check
#   #   Y_P: Y after being updated by the pilot information
#   #
#   # Note: Not actually used right now in the main function, so I'm commenting it out
#   
#   params_updated <- extract(fit_updated)
#   Y_updated <- Y_P
#   for (i in 1:length(Y_updated$mean)){
#     Y_updated$mean[i] <- mean(params_updated$theta[,i])
#   }
#   new_rank <- order(Y_updated$mean, decreasing=TRUE)
#   return(new_rank)
# }

doPilotsChangeOurMinds <- function(K,Y,original_rank,num_final_cities,Q=1){
  # Returns true if the pilots run in K change the final actions
  # and returns false otherwise, given original ranking and data Y.
  # 
  # Arguments:
  #   K: number of pilots (additional samples) run
  #   Y: set of datapoints with (mean, sd) per datpoint
  #   original_rank: int list ordering the Y datapoints from best (highest mean)
  #     to worst (lowest mean) with 1 being best
  #   num_final_cities: number of datapoints (cities) chosen in final program
  #   Q: factor by which pilot variance is different than original Y variance
  #   
  # Returns:
  #   True if the top num_final_cities are different in the original_rank vs
  #   the new rank post K cities getting additional pilots, false otherwise
    
  
  # Simulate pilots
  Y_P <- getPilotResults(K,Y,Q)
  
  # Update thetas
  updated_data <- list(I=length(Y$mean), Y=Y_P$mean, sigmaSq=Y_P$var)
  Y_updated <- extractFit(updated_data)$Y
  
  # Use new thetas to get a new city ranking
  new_rank <- order(Y_updated$mean, decreasing=TRUE)
  
  # Take the top F cities from each ranking as final city choice
  original_choice <- original_rank[1:num_final_cities]
  new_choice <- new_rank[1:num_final_cities]
  
  # Return the setequality of the two rankings (set bc order doesn't matter)
  return(!setequal(original_choice,new_choice))
}

overall<-function(data, num_pilots, num_final_cities, num_draws, Q=1){
  # Calculates the frequency at which every possible pilot set changes our minds.
  # 
  # Arguments:
  #   data: list(I,N,Y_mean,Y_sd,groups) likely generated from parseGeneratedGroupedData
  #   num_pilots: number of cities we
  #   num_final_cities: number of datapoints (cities) chosen in final program (same as K)
  #   num_draws: number of times we simulate changing our minds TODO(dsinc): better name for this?
  #   Q: factor by which pilot variance is different than original Y variance
  #   
  # Returns:
  #   True if the top num_final_cities are different in the original_rank vs
  #   the new rank post K cities getting additional pilots, false otherwise
  
  # Do initial REM fit
  Y <- extractFit(data)$Y
  
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
      nmc[i] <- nmc[i] + doPilotsChangeOurMinds(K,Y,original_rank,num_final_cities)
    }
  } 
  
  return(list(nmc=nmc, combinations=combinations))
}

data <- parseGeneratedGroupedData(generateGroupedData(I=10,mu=0,sd=10,N=3, SF=1))
# fit <- stan(file = 'randomEffectsModel2D.stan', 
#             data = data, 
#             iter = 1000, chains = 2)
# 
# extractFit(data)
# overall(data,3,2,4,Q=1)
