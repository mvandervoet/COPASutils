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
#' readSorterData("SortTest.txt", 60, 2000, 60, 5000, SVM=FALSE)
#' readSorterData("SortTest.txt", tofmin=60, extmin=60, SVM=TRUE, normalize=TRUE)

readPlate <- function(file, tofmin=0, tofmax=10000, extmin=0, extmax=10000, SVM=TRUE, normalize=FALSE) {
    plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
    modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), sort=Status.sort, TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
    modplate <- ddply(.data=modplate, .variables=c("row", "col"), .fun=extractTime)
    if(normalize){
        modplate[,10:13] <- apply(modplate[,c(5, 7:9)], 2, function(x){x/modplate$TOF})
        colnames(modplate)[10:13] <- c("norm.EXT", "norm.green", "norm.yellow", "norm.red")
    }
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
    plate = fillWells(plate)
    plot = ggplot(plate)
    if(type == "heat"){
        plot = plot + geom_rect(aes_string(xmin=0,xmax=5,ymin=0,ymax=5,fill=trait))+
               geom_text(aes_string(x=2.5,y=2.5,label=trait))+presentation+
               theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank()) +
               xlab("columns") + ylab("rows")
    }
    if(type == "scatter"){
        plot = plot + geom_point(aes_string(x = trait, y = trait2)) + presentation + xlab(trait) + ylab(trait2)
    }
    if(type == "hist"){
        plot = plot + geom_histogram(aes_string(x = trait), binwidth = diff(range(plate[[trait]]))/30) + presentation + xlab("columns") + ylab("rows")
    }
    plot = plot + facet_grid(row~col)
    return(plot)
}


#' Fill in any missing well rows with NA values
#' 
#' Returns a data frame with any missing wells filled in as NA for all measured data
#' @param plate the data frame of the plate to filled
#' @param wells the number of wells of the plate, 96 and 48 are only valid entries, defaults to 96
#' @export
#' @examples
#' fillWells(plate) #returns the plate data frame with all missing wells filled with NAs

fillWells = function(plate, wells=96){
    if(wells == 96){
        complete = data.frame(row=rep(LETTERS[1:8], each=12), col=rep(1:12, 8))
    } else if(wells == 48){
        complete = data.frame(row=rep(LETTERS[1:6], each=8), col=rep(1:8, 6))
    } else {
        stop("Invalid number of wells entered. Only 96 and 48 well plates can be filled.")
    }
    plate = merge(plate, complete, by=c("row", "col"), all=TRUE)
    return(plate)
}


#' Plot a correlation matrix within a plate or between plates
#' 
#' Returns a ggplot2 object of a correlation plot for all traits within a single plate or between two plates
#' @param plate1 one plate to compare in the correlation matrix
#' @param plate2 an optional plate to compare plate1 to, defaults to plate1, if no argument is entered plate1 will be compared to itself
#' @export
#' @examples
#' plotCorMatrix(plateA, plateB) #will plot a correlation matrix for all traits between plateA and plateB
#' plotCorMatrix(plateA) #will plot a correlation matrix for all traits within plateA

plotCorMatrix = function(plate1, plate2=plate1){
    plate1 = fillWells(plate1)
    plate2 = fillWells(plate2)
    if(nrow(plate1) != 96 | nrow(plate2) != 96){
        stop("Both plates must be summarized")
    }
    if(nrow(plate1) == nrow(plate2)){
        stop("Both plates to have the same number of traits")
    }
    corDF = melt(cor(plate1[,-(1:2)], plate2[,-(1:2)], use = "complete.obs"))
    colnames(corDF) = c("Plate1", "Plate2", "Correlation")
    ggplot(corDF, aes(Plate1, Plate2, fill = Correlation)) + 
        geom_tile() + scale_fill_gradient2("Correlation",low = "blue", high = "red", mid = "green", limits = c(-1,1)) +
        theme(axis.text.x=element_text(angle = 90, hjust = 1)) + xlab("Plate 1") + ylab("Plate 2")
}


#' Detect edge effects on 96-well plates
#' 
#' Test for an effect of the position of wells in a 96 well plate. This function will split a plate population by edge wells and center wells and test the two populations for significant differences in either a specific trait or all traits if a trait is not specified.
#' @param plate a summarized and filled plate data frame
#' @param trait a singular trait to test, defaults to NULL and will test all traits
#' @export
#' @examples
#' edgeEffect(plateA) #will return a dataframe of all pvalues for all trait tests
#' edgeEffect(plateA, "n") #will return a single value for the p value of the test with respect to the number of worms per well

edgeEffect = function(plate, trait=NULL){
    if(nrow(plate) != 96){
        stop("Plate must be summarized and filled first")
    }
    edgeWells = plate[plate$row == "A" | plate$row == "H" | plate$col == 1 | plate$col == 12,-(1:2)]
    edgeWells$pos = "edge"
    centerWells = plate[!(plate$row == "A" | plate$row == "H" | plate$col == 1 | plate$col == 12),-(1:2)]
    centerWells$pos = "center"
    total = rbind(edgeWells, centerWells)
    pos = total$pos
    total = total[,-(ncol(total))]
    if(missing(trait)){
        pval = as.data.frame(apply(total, 2, function(x){
            t.test(x[which(pos == "edge")], x[which(pos == "center")])$p.value
        }))
        pval$Trait = rownames(pval)
        rownames(pval) = NULL
        colnames(pval) = c("PValue", "Trait")
        pval = pval[,c(2,1)]
    } else {
        pval = t.test(total[,trait][which(pos == "edge")], total[,trait][which(pos == "center")])$p.value
    }
    return(pval)
}


#' Visualize and compare values and distributions across multiple plates
#' 
#' Plot the value (bar plot, if summerized) or distribution (boxplot, if unsummarized) of the data from each well across multiple plates.
#' @param plates the list of plate data frames to compare
#' @param trait the trait to compare, entered as a string
#' @param plateNames an optional character vector with the names of the individual plates; if no names are entered, numbers will be used in the order the data frames are entered
#' @export
#' @examples
#' plotCompare(list(plateA, plateB, plateC), "n", c("Plate A", "Plate B", "Plate C"))

plotCompare = function(plates, trait, plateNames=NULL){
    if(!is.null(plateNames) && length(plateNames) != length(plates)){
        stop("Length of plateNames must match length of plates")
    }
    for(i in seq(1, length(plates))){
        if(is.null(plateNames)){
            plates[[i]] = as.data.frame(cbind(plates[[i]], rep(i, nrow(plates[[i]]))))
        } else {
            plates[[i]] = as.data.frame(cbind(plates[[i]], rep(plateNames[i], nrow(plates[[i]]))))
        }
        colnames(plates[[i]])[ncol(plates[[i]])] = "Plate"
    }
    wholeDF = ldply(plates, data.frame)

    if(nrow(wholeDF) %% 96 == 0 && trait %in% colnames(wholeDF)){
        ggplot(wholeDF, aes_string(x = "Plate", y = trait, fill = "Plate")) + geom_bar(stat="identity") + facet_grid(row~col) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    } else if(trait %in% colnames(wholeDF)){
        ggplot(wholeDF, aes_string(x = "Plate", y = trait, fill = "Plate")) + geom_boxplot() + facet_grid(row~col) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    }
}
