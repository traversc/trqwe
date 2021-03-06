% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/my_functions.r
\name{trqwe_KNN}
\alias{trqwe_KNN}
\title{Fast binary KNN classifier}
\usage{
trqwe_KNN(distmat, train_idx, test_idx, classes, K, mc.cores = 1)
}
\arguments{
\item{distmat}{A NxN pre-computed distance matrix}

\item{train_idx}{Train indicies}

\item{test_idx}{Test indicies}

\item{classes}{vector length N, 1 or 0}

\item{K}{Number of nearest neighbors parameter}

\item{mc.cores}{Number of threads to use}
}
\value{
A prediction vector for the test set based on the class labels of the train set.
}
\description{
Fast binary KNN classifier
}
