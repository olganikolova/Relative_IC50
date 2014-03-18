###############################################
#
# Computing IC50 to compare drug potency
# Curve-fitting Code: courtesy of Fred Commo
# Date: 2013/09/03
#
###############################################
# Read in OncoBlau data
genpath <- "/Users/olia/Projects/_Sage/Wet-lab-collaboration/CTDD_CarlaChrisetal/Silvia_Drug_Screening/"

dir.list <- lapply(list("Oncoblau screen OVA7", "Oncoblau screen BTIC-PE","Oncoblau screening OVA6"), function(dirname){
  return(paste(genpath, dirname, sep=""))
})

file.list <- lapply(dir.list, function(dirname){
  return(paste(dirname, "/Grandori_Blau.data", sep=""))
})

data.list <- lapply(file.list, function(mypath){
  d <- read.csv(mypath, sep="\t", fill=TRUE)
  return(d)
})

ic50s.list <- lapply(as.list(seq(1:length(file.list))), function(idx){
  return(getRelativeIC50(data.list[[idx]], path=paste(dir.list[[idx]], "/Relative_IC50/", sep="")))
})

lapply(as.list(seq(1:length(file.list))), function(idx){
  fname <- paste(dir.list[[idx]], "/Relative_IC50s.csv", sep="")
  d <- do.call('rbind', ic50s.list[[idx]])
  colnames(d) <- c("Relative_IC50", "Goodness of fit")
  write.csv(d, file=fname)
})

thresholds.list <- lapply(data.list, getRespectiveThreshold) 

lapply(as.list(seq(1:length(file.list))), function(idx){
  fname <- paste(dir.list[[idx]], "/Relative_Thresholds.csv", sep="")
  d <- do.call('rbind', ic50s.list[[idx]])
  th <- do.call('rbind', thresholds.list[[idx]])
  sum(rownames(d)==rownames(th)) == nrow(d)
  res <- cbind(d,th) 
  colnames(res) <- c("Relative_IC50", "Goodness of fit","Relative threshold (0.5uM)")
  write.csv(res, file=fname)
})


###############################################
# Computes relative IC50 for a given cell line
# for a list of drugs
# Reports relative IC50 and goodness of fit
# Creates plots of the fitted curce (.pdf)
###############################################

getRelativeIC50 <- function(d=ind, 
                            path=inpath, 
                            ntpoints=8, # number of titration points
                            nreps=2     # number of replicates
                            ){
  # d <- data.list[[1]] 
  # path <- dir.list[[1]]
  # load Fred's functions
  locpath <- "/Users/olia/Projects/_Sage/Wet-lab-collaboration/CTDD_CarlaChrisetal/Silvia_Drug_Screening/Frederic/"
  source(paste(locpath,'LP5_classDef.R', sep=""))
  source(paste(locpath,'LP5_functions.R', sep=""))
  source(paste(locpath,'IC50.5P.v3_functionOnly.R', sep=""))
  
  # get uniq drug names
  drugs <- as.vector(d$Compound[d$Compound != ""])
  # fill up the first column with the corresponding drug names
  drugs.col <- rep(drugs, each=ntpoints*nreps)
  # add this to the df
  d.new <- cbind(drugs.col, d)
  # split indiv drugs
  # all replicates together
  d.list <- split(d.new, d.new$drugs.col)
  
  d.ic50s <- lapply(d.list, function(l){
     #l <-d.list$Bexarotene
    v <- l$X.Compound..M
    # standardize concentrations
    vmin <- min(v)
    vmax <- max(v)
    vrange <- vmax - vmin
    dose <- (v - vmin)/vrange
    # adjust 0s
    dose[dose==0] <- 1e-5
    dose <- log10(dose)
    
    Resp <- l$Total
    
    test <- IC50.5P(dose, Resp, Ctrl = max(Resp))
   
    # plot the graph
#     filen <- paste(path,"/", l$Compound[1], ".pdf", sep="")
#     pdf(filen)
#     plot(test)
#     dev.off()
    
    ic50 <- getEstimates(test)$D[5]
    goodnessOfFit <- summary(getGoodness(test))$adj.r.squared
    return(c(ic50,goodnessOfFit))
  
  })
  #return(d.ic50s)
} #_end_of_main_function



##########################################################
# Computes the respective concentration of relative IC50
##########################################################
getRespectiveThreshold <- function(d=ind, 
                            ntpoints=8, # number of titration points
                            nreps=2,     # number of replicates
                            thresh=0.0000005  # threshold
){
  
  # get uniq drug names
  drugs <- as.vector(d$Compound[d$Compound != ""])
  # fill up the first column with the corresponding drug names
  drugs.col <- rep(drugs, each=ntpoints*nreps)
  # add this to the df
  d.new <- cbind(drugs.col, d)
  # split indiv drugs
  # all replicates together
  d.list <- split(d.new, d.new$drugs.col)
  
  d.thresh <- lapply(d.list, function(l){
    v <- l$X.Compound..M
    # standardize concentrations
    vmin <- min(v)
    vmax <- max(v)
    vrange <- vmax - vmin
    relthresh <- (thresh - vmin)/vrange
    return(relthresh)
    
  })
  return(d.thresh)
} #_end_of_main_function
