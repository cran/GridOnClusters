//
// Created by Jiandong Wang on 2025/02/05.
//
// Copyright (c) NMSU Song lab

#include <Rcpp.h>
#include <vector>
#include "MyHash.h"
#include "evaluation.h"
#include "rcpp_utils.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List median_distance_c(const NumericMatrix & data, const NumericMatrix & labels)
{
    vector<vector<double> > c_data = convertToVector<double>(data);
    vector<vector<int> > c_labels = convertToVector<int>(labels);

    MyHash ht(c_labels);

    vector<vector<double>>  distances = evaluation_Distances(c_data, ht, "median");

    return convertToRList(distances);
}

// [[Rcpp::export]]
List mean_distance_c(const NumericMatrix & data, const NumericMatrix & labels)
{
    vector<vector<double> > c_data = convertToVector<double>(data);
    vector<vector<int> > c_labels = convertToVector<int>(labels);

    MyHash ht(c_labels);

    vector<vector<double>>  distances = evaluation_Distances(c_data, ht, "mean");

    return convertToRList(distances);
}










