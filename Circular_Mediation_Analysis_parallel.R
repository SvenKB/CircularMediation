library(foreach)
library(doParallel)
library(parallel)



#########################
#### Analyse samples ####
#########################

#### Load datasets

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
                      "/Data/Datasets",
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

simData <- loadDatasets(truen,truea,trueb,truec,nsim)

#### Analyse datasets

analyseSimDataPar <- function(simData,fun) {
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type = "FORK")
  
  
  result <- list()
  for (i in 1:length(simData)) {
    result[[i]] <- parLapply(cl,simData[[i]],fun)
  }
  estimates <- list()
  a <- matrix(0,  length(simData[[1]]),length(result[[1]]$nr1))
  for (j in 1:length(result)) {
    
    for (k in 1:length(result[[1]]$nr1)) {
      a[,k] <- as.vector(unlist(parLapply(cl,result[[j]],"[",k)))
    }
    estimates[[j]] <- a
    colnames(estimates[[j]]) <- names(result[[1]]$nr1)
  }
  stopCluster(cl)
  names(estimates) <- names(simData)
  
  return(estimates)
}


resultsDiff <- analyseSimDataPar(simData,CircMed_Diff)
resultsProd <- analyseSimDataPar(simData,CircMed_Product) 
resultsRepara <- analyseSimDataPar(simData,CircMed_Reparameter)
resultsDiffBayes <- analyseSimDataPar(simData,CircMed_Bayes_Diff)
resultsProdBayes <- analyseSimDataPar(simData,CircMed_Bayes_Product)


