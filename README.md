# run-samEst
Workshop on using the samEst R package - non-stationary spawner-recruit models for Pacific salmon

Catarina Wor (catarina.wor@dfo-mpo.gc.ca)
Dan Greenberg (dan.greenberg@dfo-mpo.gc.ca)

# Before the workshop

## Dependencies
samEst depends on the R >= 4.0.5, make sure your R version is compatible.
R > 4.0 on Windows also requires Rtools 4.X. Find out which Rtools version you need from [here](https://cran.r-project.org/bin/windows/Rtools/)

samEst depends on the R packages TMB and the latest version of stan. So we need to ensure that these packages are installed and working. To install and test the dependencies before installing samEst run [code/install_samEst.R](https://github.com/TESA-workshops/run-samEst/blob/main/code/install_samest.R)

If you run into problems, please submit an issue in this repository. 

## Exercises

Navigate to the exercises folder to find some Quarto workbooks to go through. There should be 4 exercises, the last of which will have you fitting models to real data.

## Data 
You can run the workshop on your own data, if you have it available, you will just need to format it correctly (see examples in exercises4.qmd). We are providing ~88 stocks in BC to experiment with see 'data' that you can select for exercise 4.
