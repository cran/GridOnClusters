// Created by Sajal Kumar
// Copyright (c) NMSU Song lab

#include <Rcpp.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include "Joint_Grid.h"
#include "Clusters.h"
#include "Dimension_reduction.h"
#include "evaluation.h"
#include "telescope.h"
#include "rcpp_utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List findgrid(const List cluster_info, int k, int nobs, int ndims,
                    const IntegerVector &min_bins, const IntegerVector &max_bins, String method,
                    double cutoff = 0.05, bool entropy = false, bool reduction = false, int n_thread = 1) {

    // obtain cluster centers
    // Rcpp::NumericMatrix centers = Rcpp::as<Rcpp::NumericMatrix>(cluster_info["centers"]);

    // obtain cluster labels
    Rcpp::IntegerVector labels = Rcpp::as<Rcpp::IntegerVector>(cluster_info["clusters"]);

    // obtain data
    Rcpp::NumericMatrix data = Rcpp::as<Rcpp::NumericMatrix>(cluster_info["data"]);

    // type-cast labels to std::vector<int>
    vector<int> _labels = Rcpp::as<vector<int> >(labels);

    // type-cast data to std::vector<std::vector<double > >
//    vector<vector<double> > _data(nobs, vector<double>(ndims, 0));
//    for (int i = 0; i < nobs; i++) {
//        for (int j = 0; j < ndims; j++) {
//            _data[i][j] = data(i, j);
//        }
//    }
    vector<vector<double> > _data = convertToVector<double>(data);

    // make a 'cluster' object
    //cout<<"creating cluster"<<endl;
    //sleep(5);
    Cluster clust_obj(k, _labels, _data);
    //cout<<"cluster done"<<endl;
    //sleep(5);
    
    // get grid lines
    vector<int> min_bin_limit_c = Rcpp::as<vector<int> >(min_bins);
    vector<int> max_bin_limit_c = Rcpp::as<vector<int> >(max_bins);
    //cout<<"nthread: "<<n_thread<<endl;
    Find_Grid(clust_obj, min_bin_limit_c, max_bin_limit_c, method, cutoff, entropy, reduction, n_thread);

    //cout<<"Find_Grid done"<<endl;
    //sleep(5);
    
    grid G = clust_obj.get_grids();

    Rcpp::List lines_o = convertToRList(G.lines);
    Rcpp::List table_o = convertToRList(G.conti_table);
    Rcpp::List discr_o = convertToRList(G.discr_data);
    Rcpp::List distance_to_median = convertToRList(G.median_distance);
    Rcpp::List distance_to_mean = convertToRList(G.mean_distance);
    Rcpp::NumericVector dim_bic = Rcpp::wrap(G.dimensional_bic);
    
    List lines_r = 0;
    List table_r = 0;
    List discr_r = 0;
    
    if(reduction){
        lines_r = convertToRList(G.lines_after_removed);
        table_r = convertToRList(G.conti_table_after);
        discr_r = convertToRList(G.discr_data_after);
    }

    return List::create(lines_o, table_o, discr_o, G.original_BIC, G.original_purity,
                             lines_r, table_r, discr_r, dim_bic, G.removed_BIC, G.removed_purity, 
                             G.line_removed, G.upsilon, distance_to_median, distance_to_mean);
}
