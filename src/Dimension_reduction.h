//
// Created by Jiandong Wang on 2022/11/3.
//

#ifndef GOCTEST_DIMENSION_REDUCTION_H
#define GOCTEST_DIMENSION_REDUCTION_H

#include <vector>
#include "evaluation.h"
#include "Clusters.h"
#include "discretization.h"
#include "MyHash.h"

using namespace std;

vector<bool> dim_reduct(Cluster &clusters, vector<vector<double>> &lines_removed, vector<vector<size_t>> &best_table,
                        const double global_BIC, double &best_BIC, size_t n_thread, size_t step);

double dim_BIC(Cluster &clusters, vector<vector<double>> lines_dims, vector<size_t> num_lines, size_t dim_removed,
               vector<vector<size_t>> &conti_table, size_t n_thread, size_t step);

double BIC_for_table(vector<vector<size_t>> table, size_t N, vector<size_t> nlines);

#endif //GOCTEST_DIMENSION_REDUCTION_H
