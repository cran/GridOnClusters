// Created by Jiandong Wang
// Copyright (c) NMSU Song lab

#include <vector>
#include <unistd.h>
using namespace std;

#include <Rcpp.h>
using namespace Rcpp;

#include "upsilon.h"
#include "rcpp_utils.h"

// [[Rcpp::export]]
List upsilon_c(const NumericMatrix & data, bool log) {
    vector<vector<int> > mat = convertToVector<int>(data);

    double stat;
    double df;
    double effect=0.0;

    stat = upsilon_stat(mat, df, effect);

    double pval = R::pchisq(stat, df, 0, log);

    return List::create(stat, df, pval, effect);
}

