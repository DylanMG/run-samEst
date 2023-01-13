# run-samEst
Workshop on using the samEst R package - non-stationary spawner-recruit models for Pacific salmon


# Before the workshop

## Dependencies
samEst depends on the R >= 4.0.5, make sure your R version is compatible.
R > 4.0 on Windows also requires Rtools 4.X. Find out which Rtools version you need from [here](https://cran.r-project.org/bin/windows/Rtools/)

samEst depends on the R packages TMB and the latest version of stan. So we need to ensure that these packages are installed and working. To install and test the dependencies before installing samEst run [code/install_samEst.R](https://github.com/TESA-workshops/run-samEst/blob/main/code/install_samest.R)

If you run into problems, please submit an issue in this repository. 

## Data 
You can run the workshop on your own data or examples can be provided. Let Dan Greenberg (dan.greenberg@dfo-mpo.gc.ca) and Catarina Wor (catarina.wor@dfo-mpo.gc.ca) know in advance if you would like to be assigned an example data set.


# TODO

-write scripts for the workshop
 - Read in data (your own or examples provided).
 - Run TMB models (stan would take too long??)
 - run model selection criteria 
 - plot results
 - discussion questions
