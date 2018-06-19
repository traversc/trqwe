#!/usr/bin/env bash
Rscript -e "library(Rcpp); compileAttributes('.');"
Rscript -e "library(roxygen2); roxygenise('.');"


#RRO=/usr/lib64/RRO-3.2.2/R-3.2.2/bin/R
R CMD build .
R CMD INSTALL trqwe_0.1.tar.gz
# R CMD Rd2pdf . -o "manual.pdf" --no-preview --force

git add .
git commit -am 'init'
# git push origin master
