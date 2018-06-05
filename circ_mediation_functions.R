#Circular mediation
library(circular)
library(circglmbayes)
library(coda)
library(brms)
library(rstanarm)

#Simulate data
n <- 100
a <- .5
b <- tan(.5)
c <- tan(.5)



x <- rnorm(n,0,1)
m <- rnorm(n,(a*x),1)
y <- rep(0,n)

beta <- c(c,b)
pred <- cbind(x,m)
linkfun   = function(x) 2 * atan(x)

con <- linkfun(apply(pred, 1, "%*%", beta))

y_pred <- 1+con 
err <- rvmc(n,0,5)
y <- y_pred + err
y <- as.circular(y)


plot(as.circular(y))


# Difference in coefficients

CircMed_Diff <- function(x,m,y) {
  
  # Standardize predictors
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  # Create predictor matrix
  predictors <- cbind(x,m)
  # Prepare outcome 
  outcome <- as.circular(y)
  
  # Models
  mediated_model <- circular:::lm.circular.cl(y = outcome ,x = predictors , init = c(0,0))
  direct_model <- circular:::lm.circular.cl(y = outcome,x = predictors[,1], init = 0)
  
  # Coefficients
  c <- mediated_model$coefficients[[1]]
  c_tilde <- direct_model$coefficients
  
  # Calculate effects
  direct_effect <- c
  total_effect <- c_tilde
  indirect_effect <- c_tilde - c
  
  # Effect size
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  delta <- rad2deg(total_effect) - rad2deg(direct_effect)
  
  # Prepare output
  output <- list(total_effect,direct_effect,indirect_effect,delta)
  names(output) <- c("Total","Direct","Indirect","Effect size in degrees")
  
  return(output)
  }


CircMed_Product <- function(x,m,y) {
  
  # Standardize predictors
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  # Create predictor matrix
  predictors <- cbind(x,m)
  # Prepare outcome 
  outcome <- as.circular(y)
  
  # Models
  dat <- data.frame(x,m)
  mediator_model <- lm(m~x, data=dat)
  mediated_model <- circular:::lm.circular.cl(y = outcome ,x = predictors , init = c(0,0))
  total_model <- circular:::lm.circular.cl(y = outcome,x = predictors[,1], init = 0)
  
  # Coefficients
  a <- mediator_model$coefficients[[2]]
  b <- mediated_model$coefficients[[1]]
  c <- mediated_model$coefficients[[2]]
  c_tilde <- total_model$coefficients

  # Calculate effects
  indirect_effect <- a*b
  direct_effect <- c
  total_effect <- c_tilde
  
  # Effect size
  
  
  # Prepare output
  output <- list(total_effect,direct_effect,indirect_effect)
  names(output) <- c("Total","Direct","Indirect")
  
  return(output)
  
}




CircMed_Reparameter <- function(x,m,y) {
  
  # Standardize predictors
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  # x-residualize M
  dat <- data.frame(x,m)
  resid_model <- lm(m~x, data = dat)
  m_tilde <- resid_model$residuals
  
  # Create predictor matrix
  predictors <- cbind(x,m_tilde)
  # Prepare outcome 
  outcome <- as.circular(y)
  
  # Models
  # Direct model
  mediated_model <- circular:::lm.circular.cl(y = outcome ,x = predictors , init = c(0,0))
  # a-path model
  mediator_model <- lm(m~x,data=dat)

  # Calculate effects
  total_effect <- mediated_model$coefficients[[1]]
  indirect_effect <- mediator_model$coefficients[[2]]*mediated_model$coefficients[[2]]
  direct_effect <- total_effect-indirect_effect
  
  # Prepare output
  output <- list(total_effect,direct_effect,indirect_effect)
  names(output) <- c("Total","Direct","Indirect")
  return(output)
}




CircMed_Bayes_Diff <- function(x,m,y) {
 
  # Standardize predictors
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  y <- as.circular(y)
  
  # Create dataframe
  data <- data.frame(x,m,y)
  
  
  # Models
  total_model    <- circGLM(y ~ x, data = data, Q = 50000)
  mediated_model  <- circGLM(y ~ x+m, data = data, Q = 50000)
  
  # The posterior sample of mediation effects using difference method. 
  mediation_sample_difference <- cbind(total    = total_model$bt_chain[,1], 
                                       direct   = mediated_model$bt_chain[,1],
                                       indirect = total_model$bt_chain[,1] - mediated_model$bt_chain[,1])
  
  # Obtain a summary of effects
  mcmcsum <- summary(mcmc(mediation_sample_difference)) 
  
  # Combine into a table. 
  cbind(mcmcsum$statistics, mcmcsum$quantiles)
}



CircMed_Bayes_Product <- function (x,m,y) {
  
  # Standardize predictors
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  y <- as.circular(y)
  
  # Create dataframe
  data <- data.frame(x,m,y)
  
  ## Models
  
  # Predictor-Mediator model
  #pred_med <- brm(m~x, data=data, warmup = 1000, iter = 51000, chains = 1)
  pred_med <- stan_lm(m~x, data=data, prior = NULL, chains = 2, iter = 50000, algorithm = "sampling")
  
  # Mediated model
  mediated_model  <- circGLM(y ~ x+m, data = data, Q = 50000)
  
  # Total model
  total_model <- circGLM(y ~ x, data = data, Q = 50000)

  mediation_sample_product <- cbind(total    = total_model$bt_chain[,1], 
                                       direct   = mediated_model$bt_chain[,1],
                                       indirect = mediated_model$bt_chain[,2]*posterior_samples(pred_med, pars = "x"))

  # Obtain a summary of effects
  mcmcsum <- summary(mcmc(mediation_sample_product)) 
  
  # Combine into a table. 
  cbind(mcmcsum$statistics, mcmcsum$quantiles)

  }


sim_data <- function(a,b,c,n) {
  
  x <- rnorm(n,0,1)
  m <- rnorm(n,(a*x),1)
  y <- rep(0,n)
  
  beta <- c(c,b)
  pred <- cbind(x,m)
  linkfun   = function(x) 2 * atan(x)
  
  con <- linkfun(apply(pred, 1, "%*%", beta))
  
  y_pred <- 1+con 
  err <- rvmc(n,0,5)
  y <- y_pred + err
  y <- as.circular(y)
  
  data.frame(x,m,y)
}




###################################################################################################
###################################################################################################
###################################################################################################
