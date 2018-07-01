#Circular mediation
library(circular)
library(circglmbayes)
library(coda)


# Difference in coefficients

CircMed_Diff <- function(dat) {
  
  # Standardize predictors
  x <- dat[,1]
  m <- dat[,2]
  y <- dat[,3]
  
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
  
  # Prepare output
  kappa <- mediated_model$kappa
  output <- list(total_effect,direct_effect,indirect_effect,kappa)
  names(output) <- c("Total","Direct","Indirect","Residual Kappa")
  
  return(output)
  }


CircMed_Product <- function(dat) {
  
  # Standardize predictors
  x <- dat[,1]
  m <- dat[,2]
  y <- dat[,3]
  
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
  
  # Prepare output
  kappa <- mediated_model$kappa
  output <- list(total_effect,direct_effect,indirect_effect,a,b,kappa)
  names(output) <- c("Total","Direct","Indirect","a","b","Residual kappa")
  
  return(output)
  
}




CircMed_Reparameter <- function(dat) {
  
  # Standardize predictors
  x <- dat[,1]
  m <- dat[,2]
  y <- dat[,3]
  
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
  kappa <- mediated_model$kappa
  output <- list(total_effect,direct_effect,indirect_effect,kappa)
  names(output) <- c("Total","Direct","Indirect","Resid. Kappa")
  return(output)
}




CircMed_Bayes_Diff <- function(dat) {
 
  # Standardize predictors
  x <- dat[,1]
  m <- dat[,2]
  y <- dat[,3]
  
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
  #list(cbind(mcmcsum$statistics, mcmcsum$quantiles),kappa=mediated_model$kp_mean)
  output <- list(mcmcsum$statistics[1,1],
                 mcmcsum$statistics[2,1],
                 mcmcsum$statistics[3,1],
                 mediated_model$kp_mean)
  names(output) <- c("Total","Direct","Indirect","Resid. Kappa")
  return(output)
}

CircMed_Bayes_Product <- function (dat) {
  
  # Standardize predictors
  x <- dat[,1]
  m <- dat[,2]
  y <- dat[,3]
  
  x <- (x-mean(x))/sqrt(var(x))
  m <- (m-mean(m))/sqrt(var(m))
  
  y <- as.circular(y)
  
  # Create dataframe
  data <- data.frame(x,m,y)
  
  ## Models
  
  # Predictor-Mediator model
  #pred_med <- brm(m~x, data=data, warmup = 1000, iter = 51000, chains = 1)
  pred_med <- lm(m~x, data=data)
  a <- rnorm(50000,pred_med$coefficients[2], summary(pred_med)$coefficients[2,2])
  # Mediated model
  mediated_model  <- circGLM(y ~ x+m, data = data, Q = 50000)
  
  # Total model
  total_model <- circGLM(y ~ x, data = data, Q = 50000)
  
  mediation_sample_product <- cbind(total    = total_model$bt_chain[,1], 
                                    direct   = mediated_model$bt_chain[,1],
                                    indirect = mediated_model$bt_chain[,2]*a,
                                    a = a,
                                    b = mediated_model$bt_chain[,2])
  
  # Obtain a summary of effects
  mcmcsum <- summary(mcmc(mediation_sample_product)) 
  
  # Combine into a table. 
  #list(cbind(mcmcsum$statistics, mcmcsum$quantiles),kappa=mediated_model$kp_mean)
  output <- list(mcmcsum$statistics[1,1],
                 mcmcsum$statistics[2,1],
                 mcmcsum$statistics[3,1],
                 mediated_model$kp_mean)
  names(output) <- c("Total","Direct","Indirect","Resid. Kappa")
  return(output)
}


simData <- function(a,b,c,n) {
 
  linkfun   = function(x) 2 * atan(x)

  x <- rnorm(n,0,1)
  m <- rnorm(n,(a*x),1)
  y <- rep(0,n)

  beta <- c(c,b)
  pred <- cbind(x,m)
 
  con <- linkfun(apply(pred, 1, "%*%", beta))
  
  y_pred <- 1+con 
  err <- rvmc(n,0,5)
  y <- y_pred + err
  y <- as.circular(y)
  
  data <- data.frame(x,m,y)
  
  return(data)
}

###################################################################################################
###################################################################################################
###################################################################################################