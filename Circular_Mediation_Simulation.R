#################################################
#### Simulation setup for Circular Mediation ####
#################################################

###########################
#### Preparing Samples ####
###########################

# Prepare true parameter
truea <- c(.1,.2,.4,0)
trueb <- c(.1,.2,.4,0)
truec <- c(.1,.2,.4,0-..)
truen <- c(30,100,200)
nsim=100


saveDatasets <- function(truen,truea,trueb,truec,nsim, seed = 140689) {
  set.seed(seed)
  
  # prepare all possible designs
  Alldesigns <- expand.grid(a=truea,b=trueb,c=truec,n=truen,stringsAsFactors=FALSE)
  existingDesigns <- 0
  
  for (i in 1:nrow(Alldesigns)) {
    design <- Alldesigns[i,]
    
    curr_a <- design[,1]
    curr_b <- design[,2]
    curr_c <- design[,3]
    curr_n <- design[,4]
    
    # Prepare saving datasets
    DirName <- paste0(getwd(),
                      "/Data/Datasets",
                      "n=", curr_n,
                      "a=", curr_a,
                      "b=", curr_b,
                      "c=", curr_c)
    dir.create(DirName, showWarnings = FALSE)
    
    # Simulate n datasets per design
    for (j in 1:nsim) {
      
      filename <- paste0(DirName,"/nr",j,".csv")
      
      if (!file.exists(filename)) {
        dat <- sim_data(a=curr_a,b=curr_b,c=curr_c,n=curr_n)
        write.table(dat, filename, sep = ",", row.names=FALSE, col.names=FALSE)
      } else {
        existingDesigns <- existingDesigns + 1
      }
      
    }
    
  }
  if (existingDesigns > 0) {
    cat("\n[Data generation: ", existingDesigns, "/", nsim*nrow(Alldesigns),
        " datasets already existed.]\n")
  }
}

saveDatasets(truen,truea,trueb,truec,nsim = 100)


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


