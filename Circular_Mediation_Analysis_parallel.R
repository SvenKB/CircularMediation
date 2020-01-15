library(foreach)
library(doParallel)
library(parallel)


#########################
#### Analyse samples ####
#########################

#### Load datasets


# Prepare true parameter
truea <- c(.1,.2,.4,0)
trueb <- c(.1,.2,.4,0)
truec <- c(.1,.2,.4,0)
truen <- c(30,100,200)
nsim=100

saveRDS()
readRDS()

loadDatasets <- function(truen,truea,trueb,truec,nsim) {

  
  #Prepare all possible designs
  Alldesigns <- expand.grid(a=truea,b=trueb,c=truec,n=truen,stringsAsFactors=FALSE)
  nonexistendDesigns <- 0
  data <- list()
  for (i in 1:nrow(Alldesigns)) {
    design <- Alldesigns[i,]
    
    curr_a <- design[,1]
    curr_b <- design[,2]
    curr_c <- design[,3]
    curr_n <- design[,4]
    
    # Prepare loading datasets
    DirName <- paste0(getwd(),
                      "/Data/Datasets_",
                      "n=", curr_n,
                      "a=", curr_a,
                      "b=", curr_b,
                      "c=", curr_c)
    DesignName <- paste0("n=", curr_n,
                         "a=", curr_a,
                         "b=", curr_b,
                         "c=", curr_c)
    dat <- list()
    # Load all datasets per design
    for (j in 1:nsim) {
      
      filename <- paste0(DirName,"/nr",j,".csv")
      filenumber <- paste0("nr",j)
      
      if (file.exists(filename)) {
        dat[[filenumber]] <- read.csv(filename)
        
      } else {
        nonexistendDesigns <- nonexistendDesigns + 1
      }
      
    }
    # Store list of datasets per design in a list
    data[[DesignName]] <- dat
  }
  if (nonexistendDesigns > 0) {
    cat("\n[Data loading: ", nonexistendDesigns, "/", nsim*nrow(Alldesigns),
        " datasets did not exist.]\n")
  }
  return(data)
}

simulatedData <- loadDatasets(truen,truea,trueb,truec,nsim)

#### Analyse datasets

analyseSimDataPar <- function(simData,method) {
  # Calculate the number of cores
  no_cores <- detectCores() - 2
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type = "PSOCK")
  clusterExport(cl, list("as.circular","circGLM","mcmc"))

  if (method == "CM_Diff") {
    fun <- CircMed_Diff
  } else if (method == "CM_Prod") {
    fun <- CircMed_Product
  } else if (method == "CM_Repara") {
    fun <- CircMed_Reparameter
  } else if (method == "CM_B_D") {
    fun <- CircMed_Bayes_Diff
  } else if (method == "CM_B_P") {
    fun <- CircMed_Bayes_Product
  }  
  
  result <- list()
  for (i in 1:length(simData)) {
    result[[i]] <- parLapply(cl,simData[[i]],fun)
  }
  estimates <- list()
  est <- list()
  if (method == "CM_B_D" | method == "CM_B_P" ) {
    for (i in 1:length(result)) {
      for (j in 1:length(result[[i]])) {
        est[[j]] <- rbind(cbind("Estimates"=result[[i]][[j]][[1]],result[[i]][[j]][[2]]),"Resid.Kappa"=c(result[[i]][[j]][[3]],rep(0,8)))
      }
      estimates[[i]] <- est
    }
    
  } else {
  #a <- matrix(0,  length(simData[[1]]),length(result[[1]]$nr1))
  for (j in 1:length(result)) {
    estimates[[j]] <- result[[j]]$tab
    colnames(estimates[[j]]) <- names(result[[1]]$nr1)
  }
}
  stopCluster(cl)
  names(estimates) <- names(simData)
  
  return(estimates)
}


resultsDiff <- analyseSimDataPar(simulatedData,"CM_Diff")
resultsProd <- analyseSimDataPar(simulatedData,"CM_Prod") 
resultsRepara <- analyseSimDataPar(simulatedData,"CM_Repara")
resultsDiffBayes <- analyseSimDataPar(simulatedData,"CM_B_D")
resultsProdBayes <- analyseSimDataPar(simulatedData,"CM_B_P") 

saveRDS(resultsDiff, "resultsDiff.RDS")
saveRDS(resultsProd, "resultsProd.RDS")
saveRDS(resultsRepara, "resultsRepara.RDS")
saveRDS(resultsDiffBayes, "resultsDiffBayes.RDS")
saveRDS(resultsProdBayes, "resultsProdBayes.RDS")


