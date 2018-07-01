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

simulatedData <- loadDatasets(truen,truea,trueb,truec,nsim)

#### Analyse datasets

analyseSimData <- function(simData,fun) {
  result <- list()
  for (i in 1:length(simData)) {
    result[[i]] <- lapply(simData[[i]],fun)
  }
  estimates <- list()
  a <- matrix(0,  length(simData[[1]]),length(result[[1]]$nr1))
  for (j in 1:length(result)) {
    
    for (k in 1:length(result[[1]]$nr1)) {
      a[,k] <- as.vector(unlist(lapply(result[[j]],"[",k)))
    }
    estimates[[j]] <- a
    colnames(estimates[[j]]) <- names(result[[1]]$nr1)
  }
  names(estimates) <- names(simData)
  
  return(estimates)
}


resultsDiff <- analyseSimData(simulatedData,CircMed_Diff)

resultsProd <- analyseSimData(simulatedData,CircMed_Product) 
resultsRepara <- analyseSimData(simulatedData,CircMed_Reparameter)
resultsDiffBayes <- analyseSimData(simulatedData,CircMed_Bayes_Diff)
resultsProdBayes <- analyseSimData(simulatedData,CircMed_Bayes_Product)