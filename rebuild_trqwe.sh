#!/usr/bin/env bash
cd ~/Desktop/trqwe_git
Rscript -e "library(Rcpp); compileAttributes('trqwe');"
Rscript -e "library(roxygen2); roxygenise('trqwe');"


#RRO=/usr/lib64/RRO-3.2.2/R-3.2.2/bin/R

R CMD build trqwe
R CMD INSTALL trqwe_0.1.tar.gz

R CMD Rd2pdf trqwe -o "manual.pdf"
#$RRO CMD INSTALL trqwe_0.1.tar.gz
