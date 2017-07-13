# trqwe
Performance oriented statistical metrics and utility functions in R

## Installation
```R
devtools::install_github("traversc/trqwe")
```

## A summary of functions
### Classification/model metrics and functions
* f1score - f1 classification metric
* fastAUC - Area under the receiver-operator curve
* fastPR - Precision recall curve
* fastROC - Receiver-operator curve
* mccscore - Matthew's correlation coefficient
* posteriorBalance - Adjust classification posterior probability due to unbalanced dataset
* cindex - Concordance index for survival data
* nelson_aelen_surv - Nelson Aelen estimator of cumulative hazard

### Statistical tests and functions
* ww_test - k-category Wald-Wolfowitz Runs Tests
* se - Calculate standard error
* matrixFactor - Returns a design matrix from a factor

### Misc. mathematical functions
* cosineDist - Calculates the cosine distance between rows in a matrix
* logit - Logit transformation function
* sigmoid - Calculates the sigmoid function
* topn - returns the top "n" values in a vector

### Multithreaded functions (*nix only)
* mcreadRDS - Multithreaded readRDS
* mcsaveRDS - Multithreaded saveRDS
* mcsplitapply - Apply loop over matrix subsets

### Utility functions
* tablec - Very fast tabulation of vectors writen in C++
* allDups - Like the duplicated function, but returns the first instance as well
* append<- - Append to a vector in-place
* prepend<- - Prepend to a vetor in-place
* chop - chop leading and trailing white space in a string
* head2 - Like head, but limits the number of columns returned
* tail2 - Like tail, but limits the number of columns returned
* mgrepl - Multiple grepl for a list of regular expressions
* %Q% - Infix operator for string concatenation
* TCGA_barcode - Extracts information from TCGA barcodes

### Package/environment management
* reload - Reload a package and set the position in the search order
* bioc - Easy installation bioconductor packages
* install - Easy installation of CRAN packages
* timePrompt - Adds a time prompt to the R command prompt
* statsCallback - Automated summary of R objects using callbacks
* varSizes - Returns a data.frame of all R objects and their sizes in memory