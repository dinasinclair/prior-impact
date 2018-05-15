# Set up workspace
rm(list = ls())
set.seed(1)

startTime <- Sys.time()
oldTime <- startTime

# Load packages
require(data.table)
library(ggplot2)
require(Rcpp)
require(RcppEigen) 
library(gridExtra)
library(dplyr)
library(splitstackshape)
library(sem)
library(plyr)
library(reshape)
library(matrixStats)
library(metafor)
library(stringr)

##### Helper Functions #####
stratified <- function(df, group, size, select = NULL, 
                       replace = FALSE, bothSets = FALSE) {
  if (is.null(select)) {
    df <- df
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(df)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) df[[x]] %in% select[[x]])
    df <- df[rowSums(temp) == length(select), ]
  }
  df.interaction <- interaction(df[group], drop = TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
      message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
              paste(n, collapse = ", "), "),\n.Names = c(",
              paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = replace), ])
  set1 <- do.call("rbind", temp)
  
  if (isTRUE(bothSets)) {
    set2 <- df[!rownames(df) %in% rownames(set1), ]
    list(SET1 = set1, SET2 = set2)
  } else {
    set1
  }
}

# code from previous bayesian analysis
suff_theta <- function(Y, mu, tau2, sigma2){
  #vectorize for the whole sample size instead of one sample each run 
  n <- length(Y)
  m <- length(tau2)
  mean <- (matrix(Y, ncol = n, nrow = m, byrow = TRUE)*tau2 + 
             matrix(sigma2, ncol = n, nrow = m, byrow = TRUE)*mu) / 
    (matrix(sigma2, ncol = n, nrow = m, byrow = TRUE) + tau2)
  var <- matrix(sigma2, ncol = n, nrow = m, byrow = TRUE)*tau2 /
    (matrix(sigma2, ncol = n, nrow = m, byrow = TRUE) + tau2)
  #return the mean and variance for all e_i samples
  #in matrix with nrow = # of e_i and ncol = #samples
  return(list(mean = mean, var = var))
}
# sufficient statistics for P(mu|tau, Y)
suff_mu <- function(Y, sigma2, tau2){
  #one time each calculate, hard to vectorize
  lambda <- sum(Y/(sigma2 + tau2))
  omegai <- sum(1/(sigma2 + tau2))
  omega <- 1/omegai
  #return lambda and omega defined in the reference
  return(list(lambda = lambda, omega = omega))
}

# posterior samples from P(tau|Y), also get all the suff-stats for P(mu|tau, Y)
post_tau2 <- function(Y, sigma2, tau2, n.sims){
  # P(beta|tau, Y)
  N <- length(tau2) #grid length
  log_p <- rep(NA, N)
  lambda <- log_p
  omega <- log_p
  # run all the points we choose on the tau_grid
  for(ii in 1:N){
    suff <- suff_mu(Y, sigma2, tau2[ii])
    lambda[ii] <- suff$lambda*suff$omega
    omega[ii] <- suff$omega
    log_p[ii] <- 0.5*log(omega[ii]) - 0.5*sum(log(sigma2+tau2[ii])) - 
      0.5*sum((Y - lambda[ii])^2/(sigma2 + tau2[ii]))
  }
  log_p <- log_p - max(log_p) 
  p <- exp(log_p)
  p <- p/sum(p)
  index <- sample(1:N, n.sims, replace = T, prob = p)
  tau2 <- tau2[index]
  mean <- lambda[index]
  omega <- (sqrt(omega))[index]
  # return all the selected tau2 and mean (beta) and omega (for beta). 
  return(list(tau2 = tau2, mean = mean, omega = omega))
}

#combine all the above functions for a joint sample of (theta, mu, tau)
#data is the raw data from the .csv selecting a certain group
#tau2 is in the grid
#Y <- data$treatmentcoefficient
post_mix_s <- function(data, Y, tau2, n_sim,onlymu = F){
  sigma2 <- (data$treatmentstandarderror)^2
  #first sample tau2
  tau2 <- post_tau2(Y, sigma2, tau2, n_sim)
  #also get mean and variance for beta
  mean <- tau2$mean
  omega <- tau2$omega
  tau2 <- tau2$tau2
  #calculate mu
  mu <- rnorm(n_sim)
  mu <- mu*omega + mean   
  if(onlymu){
    return(list(mu=mu))
  } else{
    #calculate suff-stats for theta
    theta <- suff_theta(Y, mu, tau2, sigma2)
    #sample theta
    theta <- array(rnorm(length(theta$mean), theta$mean, sqrt(theta$var)), 
                   dim = dim(theta$mean))
    #I2 calculation (reference page 161 "metafor.pdf")
    s2 <- (nrow(data)-1)*sum(1/sigma2)/((sum(1/sigma2))^2 - sum(1/sigma2^2))
    temp <- (mean(sqrt(tau2)))^2
    I2 <- 100*temp/(temp + s2)
    return(list(theta = t(theta), tau2 = tau2, mu = mu, I2 = I2))  
  }
  
}

normv <- function( n , mean , sd ){
  out <- rnorm( n*length(mean) , mean = mean , sd = sd )
  return( matrix( out , , ncol = n , byrow = FALSE ) )
}

##### Set control parameters #####
# Number of simulations for Bayesian routines
n_sim = 1e5

# Number of cities to implement ultimate program in
nImplements <- 3 

# Number of groups to bin cities into, set as 1 to not group cities
nGroups <- 6

# Number of permutations of groups to run (only relevant when nGroups > 1)
nGroupPermutations <- 200

# Number of draws of the simulated pilot studies to perform o estimate mean change in best treatment benefit
nDraws <- 100

# Name of folder to store outputs
outputFolder <- 'SimulationOutput'
dir.create(outputFolder, showWarnings = FALSE)


##### BEGIN MAIN CODE #####

# Load sample data from CSV file
mydata = read.csv("R_data.csv", header=TRUE)
mydata$dataID <- 1:nrow(mydata) # To track which individuals are sampled from the data set

test<-stratified(mydata, c("T_hat","group_par"), 5)

# Estimate "real" treatment effect for each city
# Used later to compare simulations' predictions of best cities to real best cities for program benefit
outOfSampleModels <- dlply(mydata[-test$dataID], "group_par", function(df) lm(change_hr ~ T_hat, data = df))
outOfSampleTreatmentEffects <- sapply(outOfSampleModels, function(x) x$coefficients[2])
save(outOfSampleModels, outOfSampleTreatmentEffects, file = file.path(outputFolder, 'outOfSampleModels.rdata'))

tslsmodel <- tsls(change_en ~ treat, ~ T_hat, data=test)
summary(tslsmodel)

models <- dlply(test, "group_par", function(df) 
  lm(change_hr ~ T_hat, data = df))

# Apply coef to each model and return a data frame
fit<-ldply(models, coef)

# Print the summary of each model
# l_ply(models, summary, .print = T)
l_ply(models, summary, .print = F)


# Treatmentcoeff = T_hat for each city
studies<-ldply(models, function(x) {r.sq <- summary(x)$r.squared
                            treatmentcoefficient <- summary(x)$coefficients[2]
                            treatmentstandarderror <- summary(x)$coefficients[,2]
                            data.frame(r.sq, treatmentcoefficient, treatmentstandarderror)})

# Remove every other line of studies data (garbage formatting)
toRemove<- seq(1, nrow(studies), 2)
studies <- studies[-toRemove, ]

#################################################################################
# first consider all cities in a random-effects model, to form priors
#################################################################################
# all one group for this model.

ds <- studies
ds$group_id <- 1

group <- unique(ds$group_id)
groupsize <- sapply(group, function(x) sum(ds$group_id==x))  
group <- group[groupsize>1] 
#estimate true theta and tau for each group using all data
TrueResult <- vector('list', length = length(group))
method="FB"
for(ii in 1:length(group)){
  job.name <- paste0("group_", group[ii])
  data <- ds[ds$group_id == group[ii], ]
  if(method == "FB"){
    gridmax <- 10*sd(data$treatmentcoefficient)
    n.grid <- 2000
    tau.grid <- seq(gridmax/n.grid, gridmax, length=n.grid) 
    #discrete distribution for tau
    tau2 <- tau.grid^2
    stat <- post_mix_s(data,Y = data$treatmentcoefficient,tau2, n_sim) 
    TrueResult[[ii]] <- matrix(c(group[ii], mean(stat$tau2), mean(stat$mu),rowMeans(stat$theta),rowVars(stat$theta)),nrow=1)  
  } 
  if(method=="EB"){ 
    stat <- rma(yi =data$treatmentcoefficient, sei = data$treatmentstandarderror, method = "EB", data = data)
    TrueResult[[ii]] <- matrix(c(group[ii], mean(stat$tau2), mean(stat$b),stat$yi),nrow=1) 
    
  }
  if(method=="FE"){ 
    stat <- rma(yi =data$treatmentcoefficient, sei = data$treatmentstandarderror, method = "FE", data = data)
    TrueResult[[ii]] <- matrix(c(group[ii], mean(stat$tau2), mean(stat$b),stat$yi),nrow=1) 
    
  }
  colnames(TrueResult[[ii]]) <- c("group_id","tau2","mu",paste0("theta",1:nrow(data)),paste0("var",1:nrow(data)))  
}
names(TrueResult) <- group  

priors<-ldply(TrueResult, rbind)



#prepare priors for sampling
keep1<-c("group_id","tau2","mu",".id","theta1","theta2","theta3","theta4","theta5","theta6","theta7","theta8","theta9","theta10","theta11","theta12","theta13","theta14","theta15","theta16","theta17","theta18","theta19","theta20","theta21","theta22","theta23","theta24","theta25","theta26","theta27","theta28","theta29","theta30","theta31")
priors1<-priors[keep1]
keep2<-c("group_id","tau2","mu",".id","var1","var2","var3","var4","var5","var6","var7","var8","var9","var10","var11","var12","var13","var14","var15","var16","var17","var18","var19","var20","var21","var22","var23","var24","var25","var26","var27","var28","var29","var30","var31")
priors2<-priors[keep2]

priors_list1 <- melt(priors1, id=c("group_id","tau2","mu",".id"))
priors_list2 <- melt(priors2, id=c("group_id","tau2","mu",".id"))
colnames(priors_list2)<-c("group_id","tau2","mu",".id","variable2", "value2")
priors_list <- cbind(priors_list1,priors_list2)

colnames(priors_list) <- c("group_id","tau2","mu",".id","city","theta_i","group2","tau22","mu2",".id2","city2","se_i")
priors_list$se_i<-sqrt(priors_list$se_i)
keep<-c("theta_i","se_i")
priors_list<-priors_list[keep]
priors_list <- cbind(priors_list, data$treatmentstandarderror)
colnames(priors_list)<-c("theta_i","se_i","sigma")

##### Determine theta* that is the "most beneficial", "true" treatment effect #####
# This is the set of cities a policymaker might decide to implement a program in, given no other studies (e.g. only looking at existing literature via meta-analysis)
theta_star <- sort(priors_list$theta_i)[1:nImplements]
city_star <- paste(sort(sapply(theta_star, FUN = function(x) which(priors_list$theta_i == x))), collapse=', ')
priorBenefit <- sum(theta_star)

##### Loop over different permutations of groupings #####
for (kk in 1:nGroupPermutations)
{
  ##### Loop over different draws of a pilot study #####
  informationValueDF <- data.frame(drawID =  1:nDraws,
                                   priorCity = city_star, # Set of cities with greatest benefit among priors
                                   priorBenefit = priorBenefit, # Aggregate benefit of those cities
                                   priorBenefitPrime_Single = NA, # Aggregeate benefit of those cities in posteriors (one model)
                                   updatedCity_Single = NA, # Set of cities with greatest benefit in posteriors (one model)
                                   updatedBenefit_Single = NA, # Aggregate benefit of those cities (one model)
                                   citiesDiff_Single = NA, # Difference in priors vs posteriors best-cities (one model), T/F
                                   benefitDiff_Single = NA, # Difference in estimated benefit (one model)
                                   pilotGroupRanking_Single = NA, # Sum-rank of treatment benefit in pilot cities (lower is better) (one model)
                                   numPilotInTop5_Single = NA, # Number of pilot cities in top 5 benefit ranking (one model)
                                   tau2_Overall = NA, # Treatment effect variance across all cities
                                   priorBenefitPrime_Grouping = NA, # Aggregeate benefit of those cities in posteriors (groupings model)
                                   updatedCity_Grouping = NA, # Set of cities with greatest benefit in posteriors (groupings model)
                                   updatedBenefit_Grouping = NA, # Aggregate benefit of those cities (groupings model)
                                   citiesDiff_Grouping = NA, # Difference in priors vs posteriors best-cities (groupings model), T/F
                                   benefitDiff_Grouping = NA, # Difference in estimated benefit (groupings model)
                                   tau2_pilot = NA, # Treatment variance of the selected pilot group
                                   pilotGroupRanking_Grouping = NA, # Sum-rank of treatment benefit in pilot cities (lower is better) (groupings model)
                                   numPilotInTop5_Grouping = NA, # Number of pilot cities in top 5 benefit ranking (groupings model)
                                   stringsAsFactors = F)
  

  posteriorsStorage_Grouping <- vector(mode = 'list', length = nDraws)
  posteriorsStorage_Single <- vector(mode = 'list', length = nDraws)
  
  # Randomly create groups
  thetaGroups <- sample(rep(1:nGroups, length.out = nrow(priors_list)))
  
  # Designate the cities where a pilot study will be simulated
  pilotCities <- which(thetaGroups == 1)
  
  # Set output file name
  outputFileName <- file.path(outputFolder, paste0('output_FB_en_', 
                                                   nGroups, 'Groups_Permutation', kk, 
                                                   '.rdata'))
  
  for (jj in 1:nDraws)
  {
    #one possible draw (in policymaker's mind) of what a study could find.
    #will need to make this into a vector later.
    #assuming the new study has a sigma 5x smaller than original:
    newStudySigma <- priors_list$sigma/10
    df<-normv(1, priors_list$theta_i, newStudySigma)
    
    df_test <- cbind(df, newStudySigma)
    colnames(df_test)<-c("treatmentcoefficient","treatmentstandarderror")
    df_test<-as.data.frame(df_test)
    df_test$group_id <- 1:nrow(df_test)
    
    
    #################################################################################
    # combine new (hypothetical) result w previous study result for each i using fixed-effects
    #################################################################################
    # this avoids double-counting, which would happen if you first combined new data with priors.
    
    ds <- studies[, c("treatmentcoefficient","treatmentstandarderror")]
    
    # all different groups for a fixed-effect model.
    ds$group_id <- 1:nrow(ds)
    both_studies <- rbind(ds,df_test)
    
    fixedEffectResult <- lapply(1:nrow(studies), function(ii){
      data <- both_studies[both_studies$group_id == ii, ]
      stat <- rma(yi =data$treatmentcoefficient, sei = data$treatmentstandarderror, method = "FE", data = data)
      #print(stat)
      
      output <- matrix(c(ii, mean(stat$tau2), mean(stat$b),mean(stat$se),stat$yi),nrow=1) 
      colnames(output) <- c("group_id","tau2","mu","se","y1","y2") 
      
      output
    })
    
    after_data <- ldply(fixedEffectResult, rbind)
    
    # the mu is what counts here rather than the theta, as it is a fixed-effect model.
    
    #################################################################################
    # combine back with other priors, but only substitute in for the values of some of the i
    #################################################################################
    
    ds <- studies
    
    # random-effects model, groups separated by closeness in theta_i priors
    ds$group_id <- thetaGroups
    
    # # random-effects model with every city in the same model
    # ds$group_id <- 1
    
    # Designate cities to conduct pilot study in
    # these are the cities whose fixed-effect model results will be combined into the random effects model
    ds[pilotCities,]$treatmentcoefficient <- after_data[pilotCities,]$mu
    ds[pilotCities,]$treatmentstandarderror <- after_data[pilotCities,]$se
    
    group <- unique(ds$group_id)
    groupsize <- sapply(group, function(x) sum(ds$group_id==x))  
    group <- group[groupsize>1] 
    group <- sort(group)
    
    
    ##### ALL CITIES IN ONE MODEL #####
    gridmax <- 10*sd(ds$treatmentcoefficient)
    n.grid <- 2000
    tau.grid <- seq(gridmax/n.grid, gridmax, length=n.grid)
    tau2 <- tau.grid^2
    
    tmpStat_Single <- post_mix_s(ds, Y = ds$treatmentcoefficient, tau2, n_sim)
    posteriors_Single <- matrix(c(mean(tmpStat_Single$tau2), 
                                  mean(tmpStat_Single$mu),
                                  rowMeans(tmpStat_Single$theta),
                                  rowVars(tmpStat_Single$theta)),nrow=1)
    colnames(posteriors_Single) <- c("tau2","mu",paste0("theta",ds$group_par),paste0("var",ds$group_par)) 
    
    # Determine new theta* after information from additional study is incorporated
    theta_updated_Single <- posteriors_Single[str_detect(colnames(posteriors_Single), 'theta')] # Vector of posterior theta values
    theta_star_updated_Single <- sort(theta_updated_Single)[1:nImplements] # Determine nImplements lowest values
    city_star_updated_Single <- paste(sort(sapply(theta_star_updated_Single, FUN = function(x) which(theta_updated_Single == x))), collapse=', ') # Determine cities corresponding to lowest thetas
    
    informationValueDF$updatedCity_Single[jj] <- city_star_updated_Single
    informationValueDF$citiesDiff_Single[jj] <- informationValueDF$updatedCity_Single[jj] != informationValueDF$priorCity[jj] # Correct form of equality?
    informationValueDF$updatedBenefit_Single[jj] <- sum(posteriors_Single[, paste0('theta', unlist(str_split(informationValueDF$updatedCity_Single[jj], pattern = ", ")))])
    
    # Update theta_hat of prior city*
    informationValueDF$priorBenefitPrime_Single[jj] <- sum(posteriors_Single[, paste0('theta', unlist(str_split(informationValueDF$priorCity[jj], pattern = ", ")))])
    
    # What is change in theta* value?
    # Always nonpositive b/c if prior selection of cities is still best, difference is zero, if it's worse, than the new selection must improve benefit
    # Or nonnegative depending on sign of benefit
    informationValueDF$benefitDiff_Single[jj] <- informationValueDF$updatedBenefit_Single[jj] - informationValueDF$priorBenefitPrime_Single[jj] 
    
    # Store tau2 value
    informationValueDF$tau2_Overall[jj] <- posteriors_Single[,'tau2']
    
    # Rank pilot group
    pilotRanks_Single <- match(pilotCities, order(theta_updated_Single))
    informationValueDF$pilotGroupRanking_Single[jj] <- sum(pilotRanks_Single)
    informationValueDF$numPilotInTop5_Single[jj] <- sum(pilotRanks_Single <= 5)
    
    posteriorsStorage_Single[[jj]] <- posteriors_Single
    
    
    
    
    ##### GROUPING MODELS #####
    #estimate true theta and tau for each group using all data
    posteriorNames_Grouping <- c(as.vector(outer(c('mu_g', 'tau2_g'), group, FUN=paste0)),
                        as.vector(outer(c('theta', 'var'), ds$group_par, FUN=paste0)))
    posteriors_Grouping <- rep(NA, length(posteriorNames_Grouping))
    names(posteriors_Grouping) <- posteriorNames_Grouping
    
    # Split data
    tmpData_Grouping <- lapply(1:length(group), FUN = function(ii) ds[ds$group_id == group[ii], ])
    
    # Run FB model
    tmpStat_Grouping <- lapply(1:length(group), FUN = function(ii) {
      gridmax <- 10*sd(tmpData_Grouping[[ii]]$treatmentcoefficient)
      n.grid <- 2000
      tau.grid <- seq(gridmax/n.grid, gridmax, length=n.grid) 
      tau2 <- tau.grid^2 #discrete distribution for tau
      
      post_mix_s(tmpData_Grouping[[ii]], Y = tmpData_Grouping[[ii]]$treatmentcoefficient, tau2, n_sim)
    })
    
    # Aggregate results
    tmpTrueResult_Grouping <- lapply(1:length(group), FUN = function(ii) {
      output <- matrix(c(group[ii],
                         mean(tmpStat_Grouping[[ii]]$tau2), 
                         mean(tmpStat_Grouping[[ii]]$mu),
                         rowMeans(tmpStat_Grouping[[ii]]$theta),
                         rowVars(tmpStat_Grouping[[ii]]$theta)),nrow=1)
      colnames(output) <- c("group_id","tau2","mu",
                            paste0("theta",tmpData_Grouping[[ii]]$group_par),
                            paste0("var",tmpData_Grouping[[ii]]$group_par)) 
      
      output
    })
    
    # Store results
    
    for(ii in 1:length(group)){
      posteriors_Grouping[paste0('mu_g', group[ii])] <- tmpTrueResult_Grouping[[ii]][,'mu']
      posteriors_Grouping[paste0('tau2_g', group[ii])] <- tmpTrueResult_Grouping[[ii]][,'tau2']
      posteriors_Grouping[paste0('theta', tmpData_Grouping[[ii]]$group_par)] <- tmpTrueResult_Grouping[[ii]][,paste0("theta",tmpData_Grouping[[ii]]$group_par)]
      posteriors_Grouping[paste0('var', tmpData_Grouping[[ii]]$group_par)] <- tmpTrueResult_Grouping[[ii]][,paste0("var",tmpData_Grouping[[ii]]$group_par)]
    }
    
    # names(TrueResult) <- group  
    
    # Determine new theta* after information from additional study is incorporated
    theta_updated_Grouping <- posteriors_Grouping[str_detect(names(posteriors_Grouping), 'theta')] # Vector of posterior theta values
    theta_star_updated_Grouping <- sort(theta_updated_Grouping)[1:nImplements] # Determine nImplements lowest values
    city_star_updated_Grouping <- paste(sort(sapply(theta_star_updated_Grouping, FUN = function(x) which(theta_updated_Grouping == x))), collapse=', ') # Determine cities corresponding to lowest thetas
    
    informationValueDF$updatedCity_Grouping[jj] <- city_star_updated_Grouping
    informationValueDF$citiesDiff_Grouping[jj] <- informationValueDF$updatedCity_Grouping[jj] != informationValueDF$priorCity[jj] # Correct form of equality?
    informationValueDF$updatedBenefit_Grouping[jj] <- sum(posteriors_Grouping[paste0('theta', unlist(str_split(informationValueDF$updatedCity_Grouping[jj], pattern = ", ")))])
    
    # Update theta_hat of prior city*
    informationValueDF$priorBenefitPrime_Grouping[jj] <- sum(posteriors_Grouping[paste0('theta', unlist(str_split(informationValueDF$priorCity[jj], pattern = ", ")))])
    
    # What is change in theta* value?
    # Always nonpositive b/c if prior selection of cities is still best, difference is zero, if it's worse, than the new selection must improve benefit
    # Or nonnegative if trying to maximize metric of benefit
    informationValueDF$benefitDiff_Grouping[jj] <- informationValueDF$updatedBenefit_Grouping[jj] - informationValueDF$priorBenefitPrime_Grouping[jj] 
    
    # Store tau2 of pilot group, hard-coded as group 1
    informationValueDF$tau2_pilot[jj] <-  posteriors_Grouping['tau2_g1']
    
    # Rank pilot group
    pilotRanks_Grouping <- match(pilotCities, order(theta_updated_Grouping))
    informationValueDF$pilotGroupRanking_Grouping[jj] <- sum(pilotRanks_Grouping)
    informationValueDF$numPilotInTop5_Grouping[jj] <- sum(pilotRanks_Grouping <= 5)
    
    posteriorsStorage_Grouping[[jj]] <- posteriors_Grouping
    
    newTime <- Sys.time()
    cat('Permutation ', kk, ', Draw ', jj, 'done after ', difftime(newTime, oldTime, units='secs'), ' seconds. Total duration: ', difftime(newTime, startTime, units='mins'), ' minutes.\n')
    oldTime <- newTime
  }
  
  posteriorsStorage_Single <- ldply(posteriorsStorage_Single, rbind)
  posteriorsStorage_Grouping <- ldply(posteriorsStorage_Grouping, rbind)
  
  save(informationValueDF, posteriorsStorage_Single, posteriorsStorage_Grouping, thetaGroups, pilotCities,
       file = outputFileName)
  
}


