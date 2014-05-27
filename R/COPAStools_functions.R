#Extract the plate number
plateno <- function(string)
{
    split<-strsplit(string,"_")
    split<-split[[1]]
    split<-split[1]
    
    split<-strsplit(split,"p")
    split<-split[[1]]
    plate<-split[2]
}

experimentName <- function(filePath){
    splitfp <- strsplit(filePath,"/")
    dirName <- splitfp[[1]][(length(splitfp[[1]]))]
    details <- strsplit(dirName,"_")[[1]][2]
    return(details)
}

#Extract the metadata info
info <- function(filePath, levels = 1){
    splitfp <- strsplit(filePath,"/")
    dirName <- splitfp[[1]][(length(splitfp[[1]])-levels)]
    
    date <- strsplit(dirName,"_")[[1]][1]
    
    details <- strsplit(dirName,"_")[[1]][2]
    
    experiment <- strsplit(details,"[0-9]+")[[1]][1]
    round <- strsplit(details,"(?i)[a-z]+")[[1]][2]
    assay <- strsplit(details,"[0-9]+")[[1]][2]
    
    split <- strsplit(splitfp[[1]][(length(splitfp[[1]]))],"_")[[1]]
    drug <- strsplit(split[2],"\\.")[[1]][1]
    plate <- strsplit(split[1],"p")[[1]][2]
    
    frame <- data.frame(date,experiment,round,assay,plate,drug)
    
    return(frame)
}

#Function to process the setup data
procSetup <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000) {
    
    #Read in the sorter data from the file
    plate <- readSorter(file, tofmin, tofmax, extmin)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column),
                                       sort = Status.sort, TOF=TOF, EXT=EXT,
                                       time=Time.Stamp, green=Green, yellow=Yellow,
                                       red=Red))
    
    #Process the data
    proc <- ddply(.data=modplate[modplate$row %in% c("A","B","C","D","E","F","G", "H"),], .var=c("row", "col"), .drop=F, .fun=function(x){
        c(pop = length(x$EXT), sorted = sum(x$sort==6), TOF = mean(x$TOF), EXT = mean(x$EXT),TOFmed=median(x$TOF),EXTmed=median(x$EXT)
        )})
    
    return(proc)
}



processWells <- function(modplate, strains=NULL, quantiles=FALSE, log=FALSE) {
    processed <- ddply(.data=modplate[modplate$call50=="object" | modplate$TOF == -1,], .variables=c("row", "col"),
                       .fun=function(x){c(n=length(x$TOF),
                                          n.sorted=sum(x$sort==6),
                                          meanTOF=mean(x$TOF),
                                          medianTOF=median(x$TOF),
                                          minTOF=as.numeric(quantile(x$TOF)[1]),
                                          if(quantiles){
                                              c(
                                              q05_TOF=as.numeric(quantile(x$TOF, probs=0.05)),
                                              q10_TOF=as.numeric(quantile(x$TOF, probs=0.1)[1]),
                                              q15_TOF=as.numeric(quantile(x$TOF, probs=0.15)[1]),
                                              q20_TOF=as.numeric(quantile(x$TOF, probs=0.2)[1]),
                                              q25_TOF=as.numeric(quantile(x$TOF, probs=0.25)[1]),
                                              q30_TOF=as.numeric(quantile(x$TOF, probs=0.3)[1]),
                                              q35_TOF=as.numeric(quantile(x$TOF, probs=0.35)[1]),
                                              q40_TOF=as.numeric(quantile(x$TOF, probs=0.4)[1]),
                                              q45_TOF=as.numeric(quantile(x$TOF, probs=0.45)[1]),
                                              q55_TOF=as.numeric(quantile(x$TOF, probs=0.55)[1]),
                                              q60_TOF=as.numeric(quantile(x$TOF, probs=0.6)[1]),
                                              q65_TOF=as.numeric(quantile(x$TOF, probs=0.65)[1]),
                                              q70_TOF=as.numeric(quantile(x$TOF, probs=0.70)[1]),
                                              q75_TOF=as.numeric(quantile(x$TOF, probs=0.75)[1]),
                                              q80_TOF=as.numeric(quantile(x$TOF, probs=0.8)[1]),
                                              q85_TOF=as.numeric(quantile(x$TOF, probs=0.85)[1]),
                                              q90_TOF=as.numeric(quantile(x$TOF, probs=0.90)[1]),
                                              q95_TOF=as.numeric(quantile(x$TOF, probs=0.95)[1])
                                              )
                                          },
                                          maxTOF=as.numeric(quantile(x$TOF)[5]),
                                          meanEXT=mean(x$EXT),
                                          medianEXT=median(x$EXT),
                                          minEXT=as.numeric(quantile(x$EXT)[1]),
                                          if(quantiles){
                                              c(
                                              q05_EXT=as.numeric(quantile(x$EXT, probs=0.05)),
                                              q10_EXT=as.numeric(quantile(x$EXT, probs=0.1)[1]),
                                              q15_EXT=as.numeric(quantile(x$EXT, probs=0.15)[1]),
                                              q20_EXT=as.numeric(quantile(x$EXT, probs=0.2)[1]),
                                              q25_EXT=as.numeric(quantile(x$EXT, probs=0.25)[1]),
                                              q30_EXT=as.numeric(quantile(x$EXT, probs=0.3)[1]),
                                              q35_EXT=as.numeric(quantile(x$EXT, probs=0.35)[1]),
                                              q40_EXT=as.numeric(quantile(x$EXT, probs=0.4)[1]),
                                              q45_EXT=as.numeric(quantile(x$EXT, probs=0.45)[1]),
                                              q55_EXT=as.numeric(quantile(x$EXT, probs=0.55)[1]),
                                              q60_EXT=as.numeric(quantile(x$EXT, probs=0.6)[1]),
                                              q65_EXT=as.numeric(quantile(x$EXT, probs=0.65)[1]),
                                              q70_EXT=as.numeric(quantile(x$EXT, probs=0.70)[1]),
                                              q75_EXT=as.numeric(quantile(x$EXT, probs=0.75)[1]),
                                              q80_EXT=as.numeric(quantile(x$EXT, probs=0.8)[1]),
                                              q85_EXT=as.numeric(quantile(x$EXT, probs=0.85)[1]),
                                              q90_EXT=as.numeric(quantile(x$EXT, probs=0.90)[1]),
                                              q95_EXT=as.numeric(quantile(x$EXT, probs=0.95)[1])
                                              )
                                          },
                                          maxEXT=as.numeric(quantile(x$EXT)[5]),
                                          mean.red=mean(x$red, na.rm=T),
                                          if(quantiles){
                                              c(
                                              q10_red=as.numeric(quantile(x$red, probs=0.1)[1]),
                                              q25_red=as.numeric(quantile(x$red, probs=0.25)[1])
                                              )
                                          },
                                          median.red=median(x$red, na.rm=T),
                                          if(quantiles){
                                              c(
                                              q75_red=as.numeric(quantile(x$red, probs=0.75)[1]),
                                              q90_red=as.numeric(quantile(x$red, probs=0.9)[1])
                                              )
                                          },
                                          mean.gr=mean(x$green, na.rm=T),
                                          median.gr=median(x$green, na.rm=T),
                                          mean.y=mean(x$yellow, na.rm=T),
                                          median.y=median(x$yellow, na.rm=T),
                                          mean.normred=mean(x$norm.red, na.rm=T),
                                          if(quantiles){
                                              c(
                                              q10_normred=as.numeric(quantile(x$norm.red, probs=0.1)[1]),
                                              q25_normred=as.numeric(quantile(x$norm.red, probs=0.25)[1])
                                              )
                                          },
                                          median.normred=mean(x$norm.red, na.rm=T),
                                          if(quantiles){
                                              c(
                                              q75_normred=as.numeric(quantile(x$norm.red, probs=0.75)[1]),
                                              q90_normred=as.numeric(quantile(x$norm.red, probs=0.9)[1]),
                                              )
                                          },
                                          if(log){
                                              c(
                                              log.meanEXT=log(mean(x$EXT)),
                                              log.medianEXT=log(median(x$EXT)),
                                              log.minEXT=as.numeric(log(quantile(x$EXT)[1])),
                                              if(quantiles){
                                                  c(
                                                  log.q05_EXT=as.numeric(log(quantile(x$EXT, probs=0.05))),
                                                  log.q10_EXT=as.numeric(log(quantile(x$EXT, probs=0.1)[1])),
                                                  log.q15_EXT=as.numeric(log(quantile(x$EXT, probs=0.15)[1])),
                                                  log.q20_EXT=as.numeric(log(quantile(x$EXT, probs=0.2)[1])),
                                                  log.q25_EXT=as.numeric(log(quantile(x$EXT, probs=0.25)[1])),
                                                  log.q30_EXT=as.numeric(log(quantile(x$EXT, probs=0.3)[1])),
                                                  log.q35_EXT=as.numeric(log(quantile(x$EXT, probs=0.35)[1])),
                                                  log.q40_EXT=as.numeric(log(quantile(x$EXT, probs=0.4)[1])),
                                                  log.q45_EXT=as.numeric(log(quantile(x$EXT, probs=0.45)[1])),
                                                  log.q55_EXT=as.numeric(log(quantile(x$EXT, probs=0.55)[1])),
                                                  log.q60_EXT=as.numeric(log(quantile(x$EXT, probs=0.6)[1])),
                                                  log.q65_EXT=as.numeric(log(quantile(x$EXT, probs=0.65)[1])),
                                                  log.q70_EXT=as.numeric(log(quantile(x$EXT, probs=0.70)[1])),
                                                  log.q75_EXT=as.numeric(log(quantile(x$EXT, probs=0.75)[1])),
                                                  log.q80_EXT=as.numeric(log(quantile(x$EXT, probs=0.8)[1])),
                                                  log.q85_EXT=as.numeric(log(quantile(x$EXT, probs=0.85)[1])),
                                                  log.q90_EXT=as.numeric(log(quantile(x$EXT, probs=0.90)[1])),
                                                  log.q95_EXT=as.numeric(log(quantile(x$EXT, probs=0.95)[1]))
                                                  )
                                              },
                                              log.maxEXT=as.numeric(log(quantile(x$EXT)[5])),
                                              log.mean.red=log(mean(x$red, na.rm=T)),
                                              if(quantiles){
                                                  c(
                                                  log.q10_red=as.numeric(log(quantile(x$red, probs=0.1)[1])),
                                                  log.q25_red=as.numeric(log(quantile(x$red, probs=0.25)[1]))
                                                  )
                                              },
                                              log.median.red=log(median(x$red, na.rm=T)),
                                              if(quantiles){
                                                  c(
                                                  log.q75_red=as.numeric(log(quantile(x$red, probs=0.75)[1])),
                                                  log.q90_red=as.numeric(log(quantile(x$red, probs=0.9)[1]))
                                                  )
                                              },
                                              log.mean.normred=log(mean(x$norm.red, na.rm=T)),
                                              if(quantiles){
                                                  c(
                                                  log.q10_normred=as.numeric(log(quantile(x$norm.red, probs=0.1)[1])),
                                                  log.q25_normred=as.numeric(log(quantile(x$norm.red, probs=0.25)[1]))
                                                  )
                                              },
                                              log.median.normred=log(mean(x$norm.red, na.rm=T)),
                                              if(quantiles){
                                                  c(
                                                  log.q75_normred=as.numeric(log(quantile(x$norm.red, probs=0.75)[1])),
                                                  log.q90_normred=as.numeric(log(quantile(x$norm.red, probs=0.9)[1]))
                                                  )
                                              }
                                              )
                                          }
                       )}, .drop=F)
    if(is.null(strains)){
        analysis <- processed
    } else {
        analysis <- data.frame(strain = as.character(strains), processed)
        analysis <- analysis[order(analysis$strain),]
        analysis <- analysis[order(analysis$row, analysis$col),]
        analysis[as.numeric(analysis$meanTOF)==-1 | is.na(analysis$meanTOF),4:ncol(analysis)] <- NA
    }
    return(analysis)
}

meltdf <- function(score){
    newscore<-data.frame(row=rep(score$row,each=1),col=rep(score$col,each=1),n=rep(score$n,each=1),f.L1=rep(score$f.L1,each=1),f.L2L3=rep(score$f.L2L3,each=1),f.L4=rep(score$f.L4,each=1),f.ad=rep(score$f.ad,each=1))
    newscore<-melt(newscore,id.var=c("row","col","n"))
    return(newscore)
}

possContam <- function(procDataFrame){
    strainMean <- mean(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    strainSD <- sd(procDataFrame$n[!is.na(procDataFrame$strain)], na.rm = TRUE)
    washMean <- mean(procDataFrame$n[is.na(procDataFrame$strain)], na.rm = TRUE)
    washSD <- sd(procDataFrame$n[is.na(procDataFrame$strain)], na.rm = TRUE)
    possibleContam <- c()
    for(j in seq(1,nrow(procDataFrame),)){
        if(!is.na(procDataFrame[j,"strain"] & !is.na(procDataFrame[j,"n"]))){
            if(procDataFrame[j,"n"] > strainMean + (2*strainSD)){
                row <- as.character(procDataFrame[j, "row"])
                col <- as.numeric(procDataFrame[j, "col"])
                adjacentWash <- procDataFrame[procDataFrame$row==row & procDataFrame$col==(col+1),"n"]
                if(adjacentWash > washMean + (2*washSD)){
                    possibleContam <- append(possibleContam, paste0(row, col))
                }
            }
        }
    } 
}

#Function to read in txt file from sorter
readSorter <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000)  {
    data <- read.delim(file=file, header=T, na.strings=c("n/a"), as.is=T, stringsAsFactors=F)
    data <- data[!is.na(data$TOF),]
    data <- data[,!is.na(data[1,])]
    data <- data[(data$TOF>=tofmin & data$TOF<=tofmax) | data$TOF == -1,]
    data <- data[(data$EXT>=extmin & data$EXT<=extmax) | data$EXT == -1,]
    data$Column <- as.factor(data$Column)
    data$Row <- as.factor(data$Row)
    return(data)
}


#Function to make time per well go from 0 to X as opposed to running for the entire plate
extractTime <- function(plate) {plate$time <- plate$time - min(plate$time); return(plate) }

#Function to take sorter raw dataframe and process to determine worm or bubble using SVM
readSorterData <- function(file, tofmin=60, tofmax=2000, extmin=20, extmax=5000, SVM=TRUE) {
    require(plyr)
    plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), sort=Status.sort, TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
    modplate <- ddply(.data=modplate, .variables=c("row", "col"), .fun=extractTime)
    if(SVM){
        require(kernlab)
        load("~/Dropbox/Biosort/Scripts and functions/bubbleSVMmodel_noProfiler.RData")
        plateprediction <- predict(bubbleSVMmodel_noProfiler, modplate[,3:length(modplate)], type="probabilities")
        modplate$object <- plateprediction[,"1"]
        modplate$call50 <- factor(as.numeric(modplate$object>0.5), levels=c(1,0), labels=c("object", "bubble"))
    }
    return(modplate)
}





#' Remove wells from a data frame representing a plate
#' 
#' Returns a data frame representing a plate with bad wells removed (either with phenotype variables set to NA or rows dropped entirely)
#' @param plate a processed data frame that has been run through process pheno
#' @param badWells a character vector consisting of all wells to remove
#' @param 

removeWells <- function(plate, badWells, drop=FALSE) {
    sp.bw <- str_split(badWells, "", 3)
    if(drop){
        for (i in seq(1, length(sp.bw))) {
            row <- sp.bw[[i]][2]
            col <- sp.bw[[i]][3]
            proc[which(plate$row == row & plate$col == col),-(1:3)] <- NA
        }
    } else {
        for (i in seq(1, length(sp.bw))) {
            row <- sp.bw[[i]][2]
            col <- sp.bw[[i]][3]
            plate = proc[(which(plate$row == row & plate$col == col)),]
        }
    }
    return(plate)
}