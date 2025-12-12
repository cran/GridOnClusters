//
// Created by Jiandong Wang on 2022/1/16.
//

#ifndef GOCTEST_CUTTING_CLUSTER_DP_H
#define GOCTEST_CUTTING_CLUSTER_DP_H

//#include "funchisq.h"
//#include "Utility.h"
#include <iostream>
#include <functional>
#include "StatDistributions.h"
#include "upsilon.h"
#include "Clusters.h"

using namespace std;

//vector<int>Cutting_Cluster_1D(vector<int> &X, unsigned int Q_max, unsigned int Q_min);

vector<int> Cutting_Cluster(vector<int> &X, unsigned int K, unsigned int Q_max, unsigned int Q_min, double cut_off);

vector<int> trace_back(unsigned int start_k, const vector<vector<int>> &P);

vector<vector<int>> prep_table(const vector<int> &X, const vector<int> &stop_posi, unsigned int Q, unsigned int K);

void print_table(const vector<vector<int>> &table);

vector<vector<int>> prep_index(vector<int> &X, int K);


#endif //GOCTEST_CUTTING_CLUSTER_DP_H
