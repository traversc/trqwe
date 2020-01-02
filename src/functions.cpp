#include <Rcpp.h>
#include <unordered_map>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]


// [[Rcpp::export]]
std::unordered_map<SEXP, int> tablec_string(CharacterVector x) {
    std::unordered_map<SEXP, int> tab;
    int n = x.size();
    for(int i=0;i<n;i++) {
        tab[x[i]]++;
    }
    return tab;
}

// [[Rcpp::export]]
std::unordered_map<int, int> tablec_int(IntegerVector x) {
    std::unordered_map<int, int> tab;
    int n = x.size();
    for(int i=0;i<n;i++) {
        tab[x[i]]++;
    }
    return tab;
}

// [[Rcpp::export]]
IntegerVector tablec_factor(IntegerVector x) {
    LogicalVector nas = is_na(x);
    bool has_nas = is_true(any(nas));
    CharacterVector _levels = x.attr("levels");
    int tab_size = has_nas ? _levels.size() + 1 : _levels.size();
    IntegerVector tab=rep(0, tab_size);
    int n = x.size();
    for(int i=0;i<n;i++) {
        if(nas[i]) {
            tab[tab_size-1]++;
        } else {
            tab[x[i]-1]++;
        }
    }
    tab.attr("names") = _levels;
    return tab;
}


// [[Rcpp::export]]
NumericMatrix fast_euclidean_dist(NumericMatrix x, NumericMatrix y) {
    int xrow = x.nrow();
    int yrow = y.nrow();
    int ndims = x.ncol();
    NumericMatrix xy(xrow, yrow);
    for(int i=0; i<xrow; i++) {
        for(int j=0; j<yrow; j++) {
            xy(i,j) = 0;
            for(int k=0; k<ndims; k++) {
                xy(i,j) += (x(i,k) - y(j,k)) * (x(i,k) - y(j,k));
            }
        }
    }
    return xy;
}
