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
    int tab_size = is_true(any(nas)) ? max(na_omit(x))+1 : max(x);
    // std::vector<int> tab(tab_size, 0);
    IntegerVector tab=rep(0, tab_size);
    int n = x.size();
    for(int i=0;i<n;i++) {
        if(nas[i]) {
            tab[tab_size-1]++;
        } else {
            tab[x[i]-1]++;
        }
    }
    tab.attr("names") = x.attr("levels");
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
