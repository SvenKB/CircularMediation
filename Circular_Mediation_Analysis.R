#########################
#### Analyse samples ####
#########################
library(boot)
#### Load datasets

# Prepare true parameter
truea <- c(.1,.2,.4,0)
trueb <- c(.1,.2,.4,0)
truec <- c(.1,.2,.4,0)
truen <- c(30,100,200)
nsim=1


?read.csv

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
        dat[[filenumber]] <- read.csv(filename,header = FALSE)
        names(dat[[filenumber]]) <- c("x","m","y")
        
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

simulatedData <- loadDatasets(truen,truea,trueb,truec,nsim=1)

#### Analyse datasets
analyseSimData <- function(simData, method) {
  
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
  if (method == "CM_Diff") {
    for (i in 1:length(simData)) {
      result[[i]] <- lapply(simData[[i]],mediationBootstrap, fun = CircMed_Diff)
    }
  } else {
    for (i in 1:length(simData)) {
      result[[i]] <- lapply(simData[[i]],fun)
    }
    
  } 
  estimates <- list()
  if (method == "CM_B_D" | method == "CM_B_P" ) {
    for (i in 1:length(result)) {
      estimates[[i]] <-list(cbind(result[[i]]$nr1[[1]],result[[i]]$nr1[[2]]),"Resid.Kappa"=result[[i]]$nr1[[3]])
    }
    
  } else {
    a <- matrix(0,  length(simData[[1]]),length(result[[1]]$nr1))
    for (j in 1:length(result)) {
    
      for (k in 1:length(result[[1]]$nr1)) {
        a[,k] <- as.vector(unlist(lapply(result[[j]],"[",k)))
      }
      estimates[[j]] <- a
      colnames(estimates[[j]]) <- names(result[[1]]$nr1)
    }
  }  
  names(estimates) <- names(simData)
  
  return(estimates)
}


testtest <- lapply(simulatedData,mediationBootstrap, fun = CircMed_Product)


resultsDiff <- analyseSimData(simulatedData,"CM_Diff")

saveRDS(resultsDiff, "resultsDiff.RDS")
resultsProd <- analyseSimData(simulatedData,"CM_Prod")
saveRDS(resultsProd, "resultsProd.RDS")
resultsRepara <- analyseSimData(simulatedData,"CM_Repara")
saveRDS(resultsRepara, "resultsRepara.RDS")
resultsDiffBayes <- analyseSimData(simulatedData,"CM_B_D")
saveRDS(resultsDiffBayes, "resultsDiffBayes.RDS")
resultsProdBayes <- analyseSimData(simulatedData,"CM_B_P") 
saveRDS(resultsProdBayes, "resultsProdBayes.RDS")

difference <- readRDS("resultsDiff.RDS")