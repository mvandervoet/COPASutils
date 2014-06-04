Using the COPASutils package to read, process, and analyze COPAS data
========================================================

In this example, we will be using the COPAStools package to read in, process, and analyze data resulting fram a Genome Wide Association Study (GWAS) using *Caenorhabditis elegans* nematode worms and a COPAS BIOSORT large particle flow cytometer. This example assumes that the COPASutils package is installed, with all necessary dependecies on your local machine. To install the COPASutils package, you can use the command `install.packages("COPASutils")`.

We will begin by requiring the COPASutils package so that we can utilize the associated functions and example data.


```r
require(COPASutils)
```

```
## Loading required package: COPASutils
## Loading required package: plyr
## Loading required package: kernlab
## Loading required package: ggplot2
```


## Reading Data

We can now read in the plate data from one of the example data sets. In the GWAS experimental design, each plate is set up with three worms sorted to each well in every other column. We have included an example data set, called "control_setup.txt" that illustrates this step of the experiment and we will read this file in first and save it to a data frame called `setupRaw`.


```r
setupRaw <- readPlate("control_setup.txt")
```


If we look at the head of this data frame we can begin to step through the different data that are output by the COPAS machine.


```r
head(setupRaw)
```

```
##   row col sort TOF EXT time green yellow red object call50
## 1   A   1    6 341 206    0    42     57  44 1.0000 object
## 2   A   1    6 345 226   62    57     61  44 1.0000 object
## 3   A   1    2  53  29  421    10      4   5 0.9987 object
## 4   A   1    6 341 185  452    60     64  47 1.0000 object
## 5   A   1    0 487 193  858    47     47  38 1.0000 object
## 6   A   1    0  58  17  920     3      3   4 0.9989 object
```


We can see that the COPAS machine groups the output data primarily by row and column of a 96-well plate, represented by the columns `row` and `col` in the above data farme. Each row in the data frame represents the readings for a single object. Working through the columns left to right, we see that the sorter returns the sort status of the object in the `sort` column (for our purposes, we are only concerned with instances where sort = 6, as these worms were sorted to the repective wells in the target plate). We then see `TOF` which stands for "time of flight" or a measure of the legth of the object in microns. Next is `EXT` or "extinction", a measure of the optical density of the object . Following this is the `time` column which represents the relative time in milliseconds???? from the first object sorted per each in well. Using this scheme, the first object to pass through the flow cell for each well with therefore be assigned a time of 0. Next are the peak height values for each of the fluorescence channels, indicated by the `green`, `yellow`, and `red` columns. Lastly the columns `object` and `call50` represent data returned from the support vector machine that probabilistically determines whether each "object" is actually an object (cell, worm, etc.) or whether it is a bubble. This feature is useful if the experiments, like ours, requires the bubble trap hardware to be bypassed. `call50` displays "bubble" if the probability of being an object (`object`) is greater than .5, or "bubble" otherwise.

If you wuld like to remove the last two columns (i.e. read in data without the help of the SVM), you can set the `SVM` argument to FALSE, as below:


```r
setupRaw2 <- readPlate("control_setup.txt", SVM = FALSE)
head(setupRaw2)
```

```
##   row col sort TOF EXT time green yellow red
## 1   A   1    6 341 206    0    42     57  44
## 2   A   1    6 345 226   62    57     61  44
## 3   A   1    2  53  29  421    10      4   5
## 4   A   1    6 341 185  452    60     64  47
## 5   A   1    0 487 193  858    47     47  38
## 6   A   1    0  58  17  920     3      3   4
```


We can also set cutoffs for minimum and maximum time of flight and extinction values as such:


```r
setupRaw3 <- readPlate("control_setup.txt", tofmin = 60, tofmax = 1000, extmin = 50, 
    extmax = 500)
head(setupRaw3)
```

```
##   row col sort TOF EXT time green yellow red object call50
## 1   A   1    6 341 206    0    42     57  44      1 object
## 2   A   1    6 345 226   62    57     61  44      1 object
## 3   A   1    6 341 185  452    60     64  47      1 object
## 4   A   1    0 487 193  858    47     47  38      1 object
## 5   A   1    0 257 248 2012    66     77  60      1 object
## 6   A   1    0 469 259 2199    63     70  57      1 object
```


Lastly, we can normalize the extinction value as well as all of the fluoresence channel values by time of flight by setting the `normalize` argument to TRUE:


```r
setupRaw4 <- readPlate("control_setup.txt", normalize = TRUE)
head(setupRaw4)
```

```
##   row col sort TOF EXT time green yellow red norm.EXT norm.green
## 1   A   1    6 341 206    0    42     57  44   0.6041    0.12317
## 2   A   1    6 345 226   62    57     61  44   0.6551    0.16522
## 3   A   1    2  53  29  421    10      4   5   0.5472    0.18868
## 4   A   1    6 341 185  452    60     64  47   0.5425    0.17595
## 5   A   1    0 487 193  858    47     47  38   0.3963    0.09651
## 6   A   1    0  58  17  920     3      3   4   0.2931    0.05172
##   norm.yellow norm.red object call50
## 1     0.16716  0.12903 1.0000 object
## 2     0.17681  0.12754 1.0000 object
## 3     0.07547  0.09434 0.9987 object
## 4     0.18768  0.13783 1.0000 object
## 5     0.09651  0.07803 1.0000 object
## 6     0.05172  0.06897 0.9989 object
```


We now see that norm.EXT, norm.green, norm.yellow, and norm.red all return the peak height value for each over the time of flight measurement.

## Processing Data



