#' Read in raw sorter data
#' 
#' Reads a raw sorter file into a dataframe, removing NA values and any objects not fitting in to the min and max cut offs.
#' @param file path to sorter data file
#' @param tofmin minimum cut off for time of flight, defaults to 0
#' @param tofmax maximum cut off for time of flight, defaults to 10000
#' @param extmin minimum cut off for extinction, defaults to 0
#' @param extmax maximum cut off for extinction, defaults to 10000
#' @export
#' @examples
#' readSorter("SortTest.txt", 60, 2000, 60, 5000)
#' readSorter("SortTest.txt", tofmin=60, extmin=60)

readSorter <- function(file, tofmin=0, tofmax=10000, extmin=0, extmax=10000)  {
    data <- read.delim(file=file, header=T, na.strings=c("n/a"), as.is=T, stringsAsFactors=F)
    data <- data[!is.na(data$TOF),]
    data <- data[,!is.na(data[1,])]
    data <- data[(data$TOF>=tofmin & data$TOF<=tofmax) | data$TOF == -1,]
    data <- data[(data$EXT>=extmin & data$EXT<=extmax) | data$EXT == -1,]
    data$Column <- as.factor(data$Column)
    data$Row <- as.factor(data$Row)
    return(data)
}


#' Set time to relative
#' 
#' Sets time relative to first well run, used in other functions, not meant to be used on its own.
#' @param plate a plate to extract the time from

extractTime <- function(plate) {plate$time <- plate$time - min(plate$time); return(plate) }


#' Read in the sorter data with minimal processing
#' 
#' Returns a minimally processed data frame, optionally with SVM-mediated bubble prediction.
#' @param file path to sorter data file
#' @param tofmin minimum cut off for time of flight, defaults to 0
#' @param tofmax maximum cut off for time of flight, defaults to 10000
#' @param extmin minimum cut off for extinction, defaults to 0
#' @param extmax maximum cut off for extinction, defaults to 10000
#' @param SVM logical dictating whether to predict bubbles with the SVM
#' @export
#' @examples
#' readSorterData("SortTest.txt", 60, 2000, 60, 5000, FALSE)
#' readSorterData("SortTest.txt", tofmin=60, extmin=60, TRUE)

readPlate <- function(file, tofmin=0, tofmax=10000, extmin=0, extmax=10000, SVM=TRUE) {
    plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), sort=Status.sort, TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
    modplate <- ddply(.data=modplate, .variables=c("row", "col"), .fun=extractTime)
    if(SVM){
        plateprediction <- predict(bubbleSVMmodel_noProfiler, modplate[,3:length(modplate)], type="probabilities")
        modplate$object <- plateprediction[,"1"]
        modplate$call50 <- factor(as.numeric(modplate$object>0.5), levels=c(1,0), labels=c("object", "bubble"))
    }
    return(modplate)
}


#' Condense all objects examined to appropriate well
#' 
#' Returns a data frame with one row per well with summary statistics.
#' @param plate plate to summarize, must be run through readSorterData, not readSorter
#' @param strains a vector of all strain (sample) names input row-wise to add to the data frame
#' @param quantiles if TRUE, columns of trait quantiles (every fifth) will be added to the output dataframe, defaults to FALSE
#' @param log if TRUE, columns of log transformed EXT and red fluorescence will be added to output dataframe, defaults to FALSE
#' @export
#' @examples
#' processWells(plate, quantiles=TRUE)
#' processWells(plate, quantiles=TRUE, log=TRUE)

summarizePlate <- function(plate, strains=NULL, quantiles=FALSE, log=FALSE) {
    require(plyr)
    processed <- suppressWarnings(ddply(.data=plate[plate$call50=="object" | plate$TOF == -1,], .variables=c("row", "col"),
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
                       )}, .drop=F))
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


#' Remove wells from a data frame representing a plate
#' 
#' Returns a data frame representing a plate with bad wells removed (either with phenotype variables set to NA or rows dropped entirely).
#' @param plate a processed data frame that has been run through processWells
#' @param badWells a character vector consisting of all wells to remove
#' @param drop a logical value dictating whether to drop wells from the data frame or set measured values to NA, defaults to FALSE
#' @export
#' @examples
#' removeWells(processedPlate, c("A1", "B7", "F12")) #Sets phenotype values for wells A1, B7, and F12 to NA
#' removeWells(processedPlate, c("A1", "B7", "F12"), drop=TRUE) #Removes rows corresponding to wells A1, B7, and F12

removeWells <- function(plate, badWells, drop=FALSE) {
    sp.bw <- str_split(badWells, "", 3)
    if(!drop){
        for (i in seq(1, length(sp.bw))) {
            row <- sp.bw[[i]][2]
            col <- sp.bw[[i]][3]
            plate[which(plate$row == row & plate$col == col),-(1:2)] <- NA
        }
    } else {
        for (i in seq(1, length(sp.bw))) {
            row <- sp.bw[[i]][2]
            col <- sp.bw[[i]][3]
            plate = plate[plate$row != row & plate$col != col,]
        }
    }
    return(plate)
}





#' Create faceted plots for every well in a 96-well plate
#' 
#' Returns ggplot2 object that is facted by row and column. By default, it will plot a heat map for the trait specified as a string. Other options include scatterplots and histograms.
#' @param plate a plate data frame, either summarized or unsummarized, to plot
#' @param trait the trait to plot in a heatmap or histogram or the independent variable in a scatterplot, enter as a string
#' @param trait2 the trait which will be the dependent variable for the scatter plot, enter as a string
#' @param type the type of plot, either "heat" for heatmap, "scatter" for scatter plot, or "hist" for histogram, defaults to "heat"
#' @export
#' @examples
#' plotTrait(plate, "n") #will create a heatmap of the number of objects per well from summarized data
#' plateTrait(plate, "TOF", "EXT", "scatter") #will create a scatter plot of EXT by TOF for each well
#' plateTrait(plate, "TOF", type="hist") #will create a histogram of TOF for each well

plotTrait = function(plate, trait, trait2=NULL, type="heat"){
    plot = ggplot(plate)
    if(type == "heat"){
        plot = plot + geom_rect(aes_string(xmin=0,xmax=5,ymin=0,ymax=5,fill=trait))+
             geom_text(aes_string(x=2.5,y=2.5,label=trait))+presentation+
             theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
    }
    if(type == "scatter"){
        plot = plot + geom_point(aes_string(x = trait, y = trait2)) + presentation
    }
    if(type == "hist"){
        plot = plot + geom_histogram(aes_string(x = trait)) + presentation
    }
    plot = plot + xlab("columns")+ylab("rows") + facet_grid(row~col)
    return(plot)
}



