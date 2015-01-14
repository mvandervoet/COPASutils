# To satisfy R CMD check --as-cran
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "Plate1", "Plate2",
    "Correlation", "strain", "dose", "variable", "value", "label", "TOF", "EXT",
    "red", "green", "yellow", "norm.EXT", "norm.red", "norm.green", "norm.yellow",
    "presentation", "bubbleSVMmodel_noProfiler"))

#' @docType package
#' @importFrom reshape2 melt
#' @importFrom stringr str_split
#' @import dplyr
#' @import knitr
#' @import kernlab
#' @import ggplot2
NULL

#' BioSorter
#' 
#' Raw data resulting from a BioSorter Machine with LP Sampler
#' @name BioSorter
NULL

#' bubbleSVMmodel_noProfiler
#' 
#' Data used by the SVM when removing bubbles in readPlate
#' @name bubbleSVMmodel_noProfiler
NULL

#' doseData
#' 
#' Data resulting from a dose response curve experiment
#' @name doseData
NULL

#' plateData1
#' 
#' Data resulting from a C. elegans GWAS experiment setup
#' @name plateData1
NULL

#' plateData2
#' 
#' Data resulting from a C. elegans GWAS experiment score
#' @name plateData2
NULL

#' presentation
#' 
#' Data for making pretty plots
#' @name presentation
NULL

#' Read in raw sorter data
#' 
#' Reads a raw sorter file into a data frame, removing NA values and any objects not fitting in to the min and max cut offs.
#' @param file path to sorter data file
#' @param tofmin minimum cut off for time of flight, defaults to 0
#' @param tofmax maximum cut off for time of flight, defaults to 10000
#' @param extmin minimum cut off for extinction, defaults to 0
#' @param extmax maximum cut off for extinction, defaults to 10000
#' @param reflx logical indicating whether ReFLx module was used (TRUE) or LP Sampler was used (FALSE), defaults to TRUE
#' @export

readSorter <- function(file, tofmin=0, tofmax=10000, extmin=0, extmax=10000, reflx=TRUE)  {
    data <- read.delim(file=file, header=T, na.strings=c("n/a"), as.is=T, stringsAsFactors=F)
    if(is.null(data$EXT) & reflx){
        stop("This file appears to have come from a machine with an LP Sampler, not a ReFLx module. Please set `reflx` = FALSE and try again.")
    }
    data <- data[!is.na(data$TOF),]
    data <- data[,!is.na(data[1,])]
    data <- data[(data$TOF>=tofmin & data$TOF<=tofmax) | data$TOF == -1,]
    if(reflx){
        data <- data[(data$EXT>=extmin & data$EXT<=extmax) | data$EXT == -1,]
        data$Column <- as.factor(data$Column)
        data$Row <- as.factor(data$Row)
    } else {
        data <- data[(data$Extinction>=extmin & data$Extinction<=extmax) | data$Extinction == -1,]
        data$Column <- as.factor(data$Column)
        data$Row <- as.factor(data$Row)
    }
    levels(data$Row) <- LETTERS[1:8]
    levels(data$Column) <- 1:12
    return(data)
}


#' Set time to relative
#' 
#' Sets time relative to first well run, used in other functions, not meant to be used on its own.
#' @param plate a plate to extract the time from
#' @export

extractTime <- function(plate){
    plate$time <- plate$time - min(plate$time)
    return(plate)
}


#' Read in the sorter data with minimal processing
#' 
#' Returns a minimally processed data frame, optionally with SVM-mediated bubble prediction.
#' @param file path to sorter data file
#' @param tofmin minimum cut off for time of flight, defaults to 0
#' @param tofmax maximum cut off for time of flight, defaults to 10000
#' @param extmin minimum cut off for extinction, defaults to 0
#' @param extmax maximum cut off for extinction, defaults to 10000
#' @param SVM logical dictating whether to predict bubbles with the SVM
#' @param reflx logical indicating whether ReFLx module was used (TRUE) or LP Sampler was used (FALSE), defaults to TRUE
#' @export

readPlate <- function(file, tofmin=0, tofmax=10000, extmin=0, extmax=10000, SVM=TRUE, reflx=TRUE) {
    plate <- readSorter(file, tofmin, tofmax, extmin, extmax, reflx)
    if(reflx){
        modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), sort=Status.sort, TOF=TOF, EXT=EXT, time=Time.Stamp, green=Green, yellow=Yellow, red=Red))
    } else {
        modplate <- with(plate, data.frame(row=Row, col=as.factor(Column), sort=Sorted.status, TOF=TOF, EXT=Extinction, time=Time, green=Green, yellow=Yellow, red=Red))
    }
    modplate <- modplate %>% group_by(row, col) %>% do(extractTime(.)) %>% data.frame()
    modplate[,10:13] <- apply(modplate[,c(5, 7:9)], 2, function(x){x/modplate$TOF})
    colnames(modplate)[10:13] <- c("norm.EXT", "norm.green", "norm.yellow", "norm.red")
    if(SVM){
        plateprediction <- predict(bubbleSVMmodel_noProfiler, modplate[,4:9], type="probabilities")
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
#' @param ends if TRUE, columns of min and max values for all traits will be added to output dataframe, defaults to FALSE
#' @export
#' @examples
#' # exampleStrains <- rep(c("N2", NA), times=48)
#' # plate1 <- summarizePlate(plateData1, quantiles=TRUE, log=TRUE, ends=TRUE)
#' # plate2 <- summarizePlate(plateData2, strains=exampleStrains)

summarizePlate <- function(plate, strains=NULL, quantiles=FALSE, log=FALSE, ends=FALSE) {
    plate <- plate[as.character(plate$call50)=="object" | plate$TOF==-1 | is.na(as.character(plate$call50)),]
    plate <- fillWells(plate)
    plate[2:14] <- lapply(plate[2:14], as.numeric)
    processed <- suppressWarnings(plate %>% group_by(row, col) %>% summarise(n=ifelse(length(TOF[!is.na(TOF)])==0, NA, length(TOF[!is.na(TOF)])),
                                                            n.sorted=sum(sort==6),
                                                            
                                                            mean.TOF=mean(TOF, na.rm=TRUE),
                                                            min.TOF=as.numeric(quantile(TOF, na.rm=TRUE)[1]),
                                                            q10.TOF=as.numeric(quantile(TOF, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.TOF=as.numeric(quantile(TOF, probs=0.25, na.rm=TRUE)[1]),
                                                            median.TOF=median(TOF, na.rm=TRUE),
                                                            q75.TOF=as.numeric(quantile(TOF, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.TOF=as.numeric(quantile(TOF, probs=0.90, na.rm=TRUE)[1]),
                                                            max.TOF=as.numeric(quantile(TOF, na.rm=TRUE)[5]),
                                                            
                                                            mean.EXT=mean(EXT, na.rm=TRUE),
                                                            min.EXT=as.numeric(quantile(EXT, na.rm=TRUE)[1]),
                                                            q10.EXT=as.numeric(quantile(EXT, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.EXT=as.numeric(quantile(EXT, probs=0.25, na.rm=TRUE)[1]),
                                                            median.EXT=median(EXT, na.rm=TRUE),
                                                            q75.EXT=as.numeric(quantile(EXT, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.EXT=as.numeric(quantile(EXT, probs=0.90, na.rm=TRUE)[1]),
                                                            max.EXT=as.numeric(quantile(EXT, na.rm=TRUE)[5]),
                                                            
                                                            mean.red=mean(red, na.rm=TRUE),
                                                            min.red=as.numeric(quantile(red, na.rm=TRUE)[1]),
                                                            q10.red=as.numeric(quantile(red, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.red=as.numeric(quantile(red, probs=0.25, na.rm=TRUE)[1]),
                                                            median.red=median(red, na.rm=TRUE),
                                                            q75.red=as.numeric(quantile(red, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.red=as.numeric(quantile(red, probs=0.9, na.rm=TRUE)[1]),
                                                            max.red=as.numeric(quantile(red, na.rm=TRUE)[5]),
                                                            
                                                            mean.green=mean(green, na.rm=TRUE),
                                                            min.green=as.numeric(quantile(green, na.rm=TRUE)[1]),
                                                            q10.green=as.numeric(quantile(green, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.green=as.numeric(quantile(green, probs=0.25, na.rm=TRUE)[1]),
                                                            median.green=median(green, na.rm=TRUE),
                                                            q75.green=as.numeric(quantile(green, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.green=as.numeric(quantile(green, probs=0.9, na.rm=TRUE)[1]),
                                                            max.green=as.numeric(quantile(green, na.rm=TRUE)[5]),
                                                            
                                                            mean.yellow=mean(yellow, na.rm=TRUE),
                                                            min.yellow=as.numeric(quantile(yellow, na.rm=TRUE)[1]),
                                                            q10.yellow=as.numeric(quantile(yellow, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.yellow=as.numeric(quantile(yellow, probs=0.25, na.rm=TRUE)[1]),
                                                            median.yellow=median(yellow, na.rm=TRUE),
                                                            q75.yellow=as.numeric(quantile(yellow, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.yellow=as.numeric(quantile(yellow, probs=0.9, na.rm=TRUE)[1]),
                                                            max.yellow=as.numeric(quantile(yellow, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.EXT=mean(norm.EXT, na.rm=TRUE),
                                                            min.norm.EXT=as.numeric(quantile(norm.EXT, na.rm=TRUE)[1]),
                                                            q10.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.EXT=median(norm.EXT, na.rm=TRUE),
                                                            q75.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.EXT=as.numeric(quantile(norm.EXT, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.red=mean(norm.red, na.rm=TRUE),
                                                            min.norm.red=as.numeric(quantile(norm.red, na.rm=TRUE)[1]),
                                                            q10.norm.red=as.numeric(quantile(norm.red, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.red=as.numeric(quantile(norm.red, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.red=median(norm.red, na.rm=TRUE),
                                                            q75.norm.red=as.numeric(quantile(norm.red, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.red=as.numeric(quantile(norm.red, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.red=as.numeric(quantile(norm.red, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.green=mean(norm.green, na.rm=TRUE),
                                                            min.norm.green=as.numeric(quantile(norm.green, na.rm=TRUE)[1]),
                                                            q10.norm.green=as.numeric(quantile(norm.green, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.green=as.numeric(quantile(norm.green, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.green=median(norm.green, na.rm=TRUE),
                                                            q75.norm.green=as.numeric(quantile(norm.green, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.green=as.numeric(quantile(norm.green, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.green=as.numeric(quantile(norm.green, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.yellow=mean(norm.yellow, na.rm=TRUE),
                                                            min.norm.yellow=as.numeric(quantile(norm.yellow, na.rm=TRUE)[1]),
                                                            q10.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.yellow=median(norm.yellow, na.rm=TRUE),
                                                            q75.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.yellow=as.numeric(quantile(norm.yellow, na.rm=TRUE)[5]),
                                                            
                                                            mean.log.EXT=mean(log(EXT), na.rm=TRUE),
                                                            min.log.EXT=as.numeric(quantile(log(EXT), na.rm=TRUE)[1]),
                                                            q10.log.EXT=as.numeric(quantile(log(EXT), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.EXT=as.numeric(quantile(log(EXT), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.EXT=median(log(EXT), na.rm=TRUE),
                                                            q75.log.EXT=as.numeric(quantile(log(EXT), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.EXT=as.numeric(quantile(log(EXT), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.EXT=as.numeric(quantile(log(EXT), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.red=mean(log(red), na.rm=TRUE),
                                                            min.log.red=as.numeric(quantile(log(red), na.rm=TRUE)[1]),
                                                            q10.log.red=as.numeric(quantile(log(red), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.red=as.numeric(quantile(log(red), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.red=median(log(red), na.rm=TRUE),
                                                            q75.log.red=as.numeric(quantile(log(red), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.red=as.numeric(quantile(log(red), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.red=as.numeric(quantile(log(red), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.green=mean(log(green), na.rm=TRUE),
                                                            min.log.green=as.numeric(quantile(log(green), na.rm=TRUE)[1]),
                                                            q10.log.green=as.numeric(quantile(log(green), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.green=as.numeric(quantile(log(green), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.green=median(log(red), na.rm=TRUE),
                                                            q75.log.green=as.numeric(quantile(log(green), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.green=as.numeric(quantile(log(green), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.green=as.numeric(quantile(log(green), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.yellow=mean(log(yellow), na.rm=TRUE),
                                                            min.log.yellow=as.numeric(quantile(log(yellow), na.rm=TRUE)[1]),
                                                            q10.log.yellow=as.numeric(quantile(log(yellow), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.yellow=as.numeric(quantile(log(yellow), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.yellow=median(log(yellow), na.rm=TRUE),
                                                            q75.log.yellow=as.numeric(quantile(log(yellow), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.yellow=as.numeric(quantile(log(yellow), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.yellow=as.numeric(quantile(log(yellow), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.EXT=mean(log(norm.EXT), na.rm=TRUE),
                                                            min.log.norm.EXT=as.numeric(quantile(log(norm.EXT), na.rm=TRUE)[1]),
                                                            q10.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.EXT=median(log(norm.EXT), na.rm=TRUE),
                                                            q75.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.EXT=as.numeric(quantile(log(norm.EXT), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.red=mean(log(norm.red), na.rm=TRUE),
                                                            min.log.norm.red=as.numeric(quantile(log(norm.red), na.rm=TRUE)[1]),
                                                            q10.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.red=median(log(norm.red), na.rm=TRUE),
                                                            q75.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.red=as.numeric(quantile(log(norm.red), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.green=mean(log(norm.green), na.rm=TRUE),
                                                            min.log.norm.green=as.numeric(quantile(log(norm.green), na.rm=TRUE)[1]),
                                                            q10.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.green=median(log(norm.red), na.rm=TRUE),
                                                            q75.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.green=as.numeric(quantile(log(norm.green), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.yellow=mean(log(norm.yellow), na.rm=TRUE),
                                                            min.log.norm.yellow=as.numeric(quantile(log(norm.yellow), na.rm=TRUE)[1]),
                                                            q10.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.yellow=median(log(norm.yellow), na.rm=TRUE),
                                                            q75.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.yellow=as.numeric(quantile(log(norm.yellow), na.rm=TRUE)[5])))

    if(!ends){
        processed <- processed[,-(grep("min", colnames(processed)))]
        processed <- processed[,-(grep("max", colnames(processed)))]
    }
    if(!quantiles){
        processed <- processed[,-(grep("q", colnames(processed)))]
    }
    if(!log){
        processed <- processed[,-(grep("log", colnames(processed)))]
    }
    if(is.null(strains)){
        analysis <- processed
        analysis <- analysis[order(analysis$row, analysis$col),]
    } else {
        analysis <- data.frame(strain = as.character(strains), processed)
        analysis <- analysis[order(analysis$strain),]
        analysis <- analysis[order(analysis$row, analysis$col),]
    }
    analysis[analysis$mean.TOF==-1 | is.na(analysis$mean.TOF),which(colnames(analysis)=="n"):ncol(analysis)] <- NA
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
#' # plate <- summarizePlate(plateData1)
#' # plateWithoutWells <- removeWells(plate=plate, badWells=c("A1", "A2", "A3"))

removeWells <- function(plate, badWells, drop=FALSE) {
    sp.bw <- str_split(badWells, "", 3)
    if(!drop){
        for (i in seq(1, length(sp.bw))) {
            if(length(sp.bw) > 0){
                row <- as.character(sp.bw[[i]][2])
                col <- as.character(sp.bw[[i]][3])
                plate[which(plate$row == row & plate$col == col),which(colnames(plate)=="n"):ncol(plate)] <- NA
            }
        }
    } else {
        for (i in seq(1, length(sp.bw))) {
            if(length(sp.bw) > 0){
                row <- sp.bw[[i]][2]
                col <- sp.bw[[i]][3]
                plate = plate[plate$row != row | plate$col != col,]
            }
        }
    }
    return(plate)
}


#' Create faceted plots for every well in a 96-well plate
#' 
#' Returns ggplot2 object that is faceted by row and column. By default, it will plot a heat map for the trait specified as a string. Other options include scatterplots and histograms.
#' @param plate a plate data frame, either summarized or unsummarized, to plot
#' @param trait the trait to plot in a heat map or histogram or the independent variable in a scatter plot, enter as a string
#' @param trait2 the trait which will be the dependent variable for the scatter plot, enter as a string
#' @param type the type of plot, either "heat" for heatmap, "scatter" for scatter plot, or "hist" for histogram, defaults to "heat"
#' @export
#' @examples
#' #### COPASutils Figures
#' ### Figure 1a 
#' #plotTrait(doseData, trait="TOF", trait2="EXT", type="scatter")
#' 
#' ### Figure 1b
#' #plotTrait(doseData, trait="TOF", type="hist")
#' 
#' ### Figure 1c
#' # sumDose <- summarizePlate(doseData)
#' # plotTrait(sumDose, trait="n", type="heat") 

plotTrait = function(plate, trait, trait2=NULL, type="heat"){
    plate = data.frame(plate)
    if(type == "heat"){
        plate$label <- ifelse(is.na(plate[,which(colnames(plate)==trait)]), "", plate[,which(colnames(plate)==trait)])
        plot = ggplot(plate) + geom_rect(aes_string(xmin=0,xmax=5,ymin=0,ymax=5,fill=trait)) +
               geom_text(aes(x=2.5,y=2.5,label=label, colour="white"), colour="white")+presentation +
               theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank()) +
               xlab("Columns") + ylab("Rows")
    } else if(type == "scatter"){
        plot = ggplot(plate) + geom_point(aes_string(x = trait, y = trait2), size=1.5) + presentation + xlab(trait) + ylab(trait2) +
               theme(axis.text.x=element_text(size="10", angle=45, hjust=1), axis.text.y=element_text(size="10"))
    } else if(type == "hist"){
        if(sum(is.na(plate[,which(colnames(plate)==trait)])) > 0){
            badWells = plate[is.na(plate[,which(colnames(plate)==trait)]),c("row", "col")]
            badWells = paste0(badWells$row, badWells$col)
            plate = removeWells(plate, badWells, drop=TRUE)
        }
        plot = ggplot(plate) + geom_histogram(aes_string(x = trait), binwidth = diff(range(plate[[trait]]))/15) + presentation + xlab(trait) + ylab("Count") +
               theme(axis.text.x=element_text(size="10", angle=45, hjust=1), axis.text.y=element_text(size="10"))
    } else {
        stop("Unrecognized plot type")
    }
    plot = plot + facet_grid(row~col, drop=FALSE)
    return(plot)
}


#' Fill in any missing well rows with NA values
#' 
#' Returns a data frame with any missing wells filled in as NA for all measured data
#' @param plate the data frame of the plate to filled
#' @export

fillWells = function(plate){
    complete = data.frame(row=rep(LETTERS[1:8], each=12), col=rep(1:12, 8))
    plate = merge(plate, complete, by=c("row", "col"), all=TRUE)
    return(plate)
}


#' Plot a correlation matrix within a plate or between plates
#' 
#' Returns a ggplot2 object of a correlation plot for all traits within a single plate or between two plates. The data from each plate must be summarized prior to being passed to this function.
#' @param plate_summary1 one summarized plate to compare in the correlation matrix
#' @param plate_summary2 an optional summarized plate to compare with plate_summary1, defaults to plate_summary1; if no argument is entered plate_summary1 will be compared to itself
#' @export
#' @examples
#' #### COPASutils Figures
#' ### Figure 4
#' # library(dplyr)
#' # sumDose <- summarizePlate(doseData)
#' # corDose <- select(sumDose, -n.sorted)
#' # plotCorMatrix(corDose)

plotCorMatrix = function(plate_summary1, plate_summary2=plate_summary1){
    plate_summary1 = fillWells(plate_summary1)
    plate_summary2 = fillWells(plate_summary2)
    if(nrow(plate_summary1) != 96 | nrow(plate_summary2) != 96){
        stop("Both plates must be summarized")
    }
    if(ncol(plate_summary1) != ncol(plate_summary2)){
        stop("Both plates have to have the same number of traits")
    }
    corDF = melt(cor(plate_summary1[,-(1:2)], plate_summary2[,-(1:2)], use = "complete.obs"))
    colnames(corDF) = c("Plate1", "Plate2", "Correlation")
    ggplot(corDF, aes(Plate1, Plate2, fill = Correlation)) + 
        geom_tile() + scale_fill_gradient2("Correlation",low = "blue", high = "red", mid = "green", limits = c(-1,1)) +
        theme(axis.text.x=element_text(angle = 90, hjust = 1)) + xlab("Plate 1") + ylab("Plate 2")
}


#' Detect edge effects on 96-well plates
#' 
#' Test for an effect of the position of wells in a 96 well plate. This function will split a plate population by edge wells and center wells and test the two populations for significant differences in either a specific trait or all traits if a trait is not specified.
#' @param plate_summary a summarized and filled plate data frame
#' @param trait a singular trait to test, defaults to NULL and will test all traits
#' @export

edgeEffect = function(plate_summary, trait=NULL){
    plate_summary = data.frame(plate_summary)
    if(nrow(plate_summary) != 96){
        stop("plate_summary must be summarized and filled first")
    }
    edgeWells = plate_summary[plate_summary$row == "A" | plate_summary$row == "H" | plate_summary$col == 1 | plate_summary$col == 12,-(1:2)]
    edgeWells$pos = "edge"
    centerWells = plate_summary[!(plate_summary$row == "A" | plate_summary$row == "H" | plate_summary$col == 1 | plate_summary$col == 12),-(1:2)]
    centerWells$pos = "center"
    total = rbind(edgeWells, centerWells)
    pos = total$pos
    total = total[,-(ncol(total))]
    if(missing(trait)){
        pval = as.data.frame(apply(total, 2, function(x){
            wilcox.test(x[which(pos == "edge")], x[which(pos == "center")])$p.value
        }))
        pval$Trait = rownames(pval)
        rownames(pval) = NULL
        colnames(pval) = c("PValue", "Trait")
        pval = pval[,c(2,1)]
    } else {
        pval = wilcox.test(total[,trait][which(pos == "edge")], total[,trait][which(pos == "center")])$p.value
    }
    return(pval)
}


#' Visualize and compare values and distributions across multiple plates
#' 
#' Plot the value (bar plot, if summarized) or distribution (boxplot, if unsummarized) of the data from each well across multiple plates.
#' @param plates the list of plate data frames to compare
#' @param trait the trait to compare, entered as a string
#' @param plateNames an optional character vector with the names of the individual plates; if no names are entered, numbers will be used in the order the data frames are entered
#' @export
#' @examples
#' #### COPASutils Figures
#' ### Figure 3a
#' # plotCompare(list(doseData, plateData2), "TOF")
#' 
#' ### Figure 3b
#' # sumDose <- summarizePlate(doseData)
#' # sumPlate <- summarizePlate(plateData2)
#' # plotCompare(list(sumDose, sumPlate), "n")

plotCompare = function(plates, trait, plateNames=NULL){
    if(!is.null(plateNames) & length(plateNames) != length(plates)){
        stop("Length of plateNames must match length of plates")
    }
    for(i in seq(1, length(plates))){
        if(is.null(plateNames)){
            plates[[i]] = as.data.frame(cbind(plates[[i]], as.factor(rep(i, nrow(plates[[i]])))))
        } else {
            plates[[i]] = as.data.frame(cbind(plates[[i]], rep(plateNames[i], nrow(plates[[i]]))))
        }
        colnames(plates[[i]])[ncol(plates[[i]])] = "Plate"
    }
    wholeDF = suppressWarnings(rbind_all(plates))
    wholeDF[,c("row", "col")] <- list(as.factor(wholeDF$row), as.factor(wholeDF$col))
    levels(wholeDF$row) <- LETTERS[1:8]
    levels(wholeDF$col) <- 1:12

    if(nrow(wholeDF) %% 96 == 0 & trait %in% colnames(wholeDF)){
        ggplot(wholeDF, aes_string(x = "Plate", y = trait, fill = "Plate")) + geom_bar(stat="identity") + facet_grid(row~col, drop=FALSE) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    } else if(trait %in% colnames(wholeDF)){
        ggplot(wholeDF, aes_string(x = "Plate", y = trait, fill = "Plate")) + geom_boxplot() + facet_grid(row~col, drop=FALSE) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    }
}

#' Plot a dose response curve by strain across a plate
#' 
#' Return a ggplot2 object of a dose response curve by strain across a plate
#' @param plate_summary the summarized plate data frame to be plotted
#' @param dosages a vector of dosages in the plate, entered by row
#' @param trait the trait to be plotted on the y axis
#' @export
#' @examples
#' #### COPASutils figure
#' ### Figure 5
#' # strains <- rep(c("A", "B", "C", "D"), each=6, times=4)
#' # sumDose <- summarizePlate(doseData, strains=strains)
#' # dose <- rep(c(1, 5, 10, 15, 20, NA), times=16)
#' # plotDR(sumDose, "n", dosages=dose)

plotDR = function(plate_summary, dosages, trait="n"){
    plate_summary = data.frame(cbind(dose = dosages, plate_summary))
    plate_summary = plate_summary[,-c(3,4)]
    plate_summary = plate_summary[!is.na(plate_summary$strain)&!is.na(plate_summary$dose),]
    plate_summary = melt(plate_summary, id=c("strain", "dose"))
    plate_summary = plate_summary %>% group_by(strain, dose, variable) %>% summarize(mean = mean(value))
    plate_summary = reshape2::dcast(plate_summary, strain+dose~variable, value.var="mean")
    ggplot(plate_summary, aes_string(x="dose", y=trait, colour="strain")) + geom_point() + geom_line() + scale_colour_discrete(name="Strain") + xlab("Dose") + ylab(trait)
}

#' Plot dose response curves for all traits
#' 
#' Return a list of ggplot2 objects of dose response curves by strain across a plate. Plots for specific traits can be accessed by name from the returned list (i.e. "plots$n" will return the dose response plot for the trait "n").
#' @param plate_summary the summarized plate data frame to be plotted
#' @param dosages a vector of dosages in the plate, entered by row
#' @export

plotDR_allTraits = function(plate_summary, dosages){
    plots = list()
    for(trait in colnames(plate_summary)[4:ncol(plate_summary)]){
        plots = append(plots, list(plotDR(plate_summary, dosages, trait)))
    }
    names(plots) = colnames(plate_summary)[4:ncol(plate_summary)]
    return(plots)
}
