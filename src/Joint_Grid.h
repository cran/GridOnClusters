//
// Created by Jiandong Wang on 2/14/20.
//
// Copyright (c) NMSU Song lab


#ifndef JOINT_GRID_JOINT_GRID_H
#define JOINT_GRID_JOINT_GRID_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <climits>
#include <cassert>
#include "Clusters.h"
#include "evaluation.h"
#include "Cut_by_density.h"
#include "Cutting_Cluster_dp.h"
#include "Cutting_Cluster_dp_compressed.h"
#include "discretization.h"
#include "MyHash.h"
#include "Dimension_reduction.h"
#include "upsilon.h"
#include "utils.h"

using namespace std;

struct grid Find_Grid(Cluster& clusters, const vector<int>& min_bin_limit, const vector<int>& max_bin_limit,
                      const string& method = "DP approx likelihood", double cut_off = 0.05,
                      bool entropy = false, bool dim_reduct = false, size_t n_thread = 1);

void Find_Grid_thread(Cluster& clusters, const vector<int>& min_bin_limit, const vector<int>& max_bin_limit, const string& method,
                      double cut_off, bool entropy, vector<vector<double>>& lines_multi_dim, vector<size_t>& num_lines,
                      vector<vector<double>>& medians_multi_dim, size_t dim_u, size_t dim_l, vector<double> & best_bic_all);

vector<double>
Find_1D_Grid(Cluster& clusters, int dim_input, int min_bin, int max_bin, vector<double>& medians, const string& method,
double cut_off, bool entropy, double & best_bic, const bool& increasing = true);

bool comp(vector<double>& x, vector<double>& y);

//static const std::vector<double> *spx;

// static bool compi(size_t i, size_t j);

bool test_sorted(const vector<double>& x, bool increasing);

double euclidean_distance(const vector<double>& point1, const vector<double>& point2);
#endif //JOINT_GRID_JOINT_GRID_H
