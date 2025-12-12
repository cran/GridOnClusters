//
// Created by Jiandong Wang on 2022/2/23.
//

#ifndef GOCTEST_CUT_BY_DENSITY_H

#include <vector>
#include <limits>
#include "Clusters.h"

#define GOCTEST_CUT_BY_DENSITY_H


vector<double> Cut_by_density_1D(Cluster &clusters, int dim_input, int min_bin_limit);

vector<vector<double> > prep_index(vector<double> &c1, vector<double> &c2, double middle_1, double middle_2);

double binary_search_index(const vector<vector<double> > &c_index, const int left, const int right, const int size_c1,
                           const int size_c2, bool &overlap, vector<double> &err_sum);


#endif //GOCTEST_CUT_BY_DENSITY_H
