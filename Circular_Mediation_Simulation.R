# Loading needed packages
library(circular)
library(circglmbayes)
library(coda)
library(brms)
library(rstanarm)


# Setting up simulation paramter

iter <- 4

# Determine Sample Size
n <- 100

# Define effects
a <- .5
b <- tan(.5)  # using tangens, so it is easier to define the effect on circular outcome
c <- tan(.5)

# Define link function
linkfun   = function(x) 2 * atan(x)

# Prepate output matrix

prod_results <- matrix(nrow = iter,ncol=3,dimnames = list(1:iter, c("Total","Direct","Indirect")))
diff_results <- matrix(nrow = iter,ncol=3,dimnames = list(1:iter, c("Total","Direct","Indirect")))
para_results <- matrix(nrow = iter,ncol=3,dimnames = list(1:iter, c("Total","Direct","Indirect")))
bayesdiff_results <- matrix(nrow = iter,ncol=3,dimnames = list(1:iter, c("Total","Direct","Indirect")))
bayesprod_results <- matrix(nrow = iter,ncol=3,dimnames = list(1:iter, c("Total","Direct","Indirect")))

  
for (i in 1:iter) {
  
  ## Simulate data
  
  x <- rnorm(n,0,1)
  m <- rnorm(n,(a*x),1)
  
  # Prepare data to simulate circular variable
  y <- rep(0,n)
  beta <- c(c,b)
  pred <- cbind(x,m)
  
  # Using link function and matrix multiplication to obtain preducited y-value
  con <- linkfun(pred %*% as.matrix(beta))
  
  # Add some noise to the predicted y-value
  y_pred <- 1+con 
  err <- rvmc(n,0,1)
  y <- y_pred + err
  
  ## Analyse data using all 5 methods
  prod <- CircMed_Product(x,m,y)
  diff <- CircMed_Diff(x,m,y)
  para <- CircMed_Reparameter(x,m,y)
  bayesdiff <- CircMed_Bayes_Diff(x,m,y)
  bayesprod <- CircMed_Bayes_Product(x,m,y)
  
  ## Save statistics of interest ##
  
  # Save output
  prod_results[i,] <- unlist(prod)
  diff_results[i,] <- unlist(diff[1:3])
  para_results[i,] <- unlist(para)
  bayesdiff_results[i,] <- bayesdiff[,1]
  bayesprod_results[i,] <-  bayesprod[,1]

}


apply(prod_results,2,mean)
apply(diff_results,2,mean)
apply(para_results,2,mean)
apply(bayesdiff_results,2,mean)
apply(bayesprod_results,2,mean)
