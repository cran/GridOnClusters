/**
 * @file Cutting_Cluster_dp_compressed.h
 * @author Jiandong Wang
 * @date 3/08/22
 * @copyright NMSU Song lab
 */

#ifndef GOCTEST_CUTTING_CLUSTER_DP_COMPRESSED_H
#define GOCTEST_CUTTING_CLUSTER_DP_COMPRESSED_H

#include <algorithm>
#include <cfloat>
#include <limits>
#include <iomanip>
//#include "upsilon.h"
#include "Clusters.h"
#include "Cutting_Cluster_dp.h"
#include "evaluation.h"
//#include "StatDistributions.h"

using namespace std;

vector<int>
Cutting_Cluster_dp_compressed(vector<int> &labels, vector<double> &data, size_t K, size_t Q_max,
                              size_t Q_min, double & best_bic, const string &method, const bool entropy);

void compress_labels(const vector<int> &labels, vector<int> &labels_compressed, vector<int> &labels_weights,
                     vector<int> &weight_sum, vector<vector<int>> &index, const unsigned int N, const unsigned int K);

void update_index(vector<vector<int>> &index, int cur_index, int cur_weight, int cur_label, unsigned int K);

vector<vector<int>> prep_table_dp(const vector<int> &label, const vector<int> &weight, const vector<int> &stop_posi,
                                  unsigned int Q, unsigned int K);

vector<int> trace_back_dp(unsigned int start_k, const vector<vector<int>> &P);

void
cal_score_absolute(const vector<int> &labels_compressed, const vector<int> &labels_weights, vector<vector<double>> &S,
                   vector<vector<int>> &P, const vector<vector<int>> &index, unsigned int Q_max, unsigned int K,
                   int n_blocks);

void cal_score_binomial(const vector<double> &data, const vector<int> &labels_compressed,
                        const vector<int> &labels_weights, vector<vector<double>> &S, vector<vector<double>> &S_e,
                        vector<vector<int>> &P, const vector<vector<int>> &index, unsigned int Q_max, unsigned int K,
                        int n_blocks, bool entropy);

void cal_score_binomial_recur_main(const vector<double> &data, const vector<int> &labels_compressed,
                                   const vector<int> &labels_weights,
                                   vector<vector<double>> &S, vector<vector<int>> &P, const vector<vector<int>> &index,
                                   unsigned int Q_max, unsigned int K, int n_blocks, bool entropy);

void cal_score_binomial_recur(int q, int il, int jl, int iu, int ju, vector<vector<double>> &S,
                              vector<vector<int>> &J, const vector<vector<int>> &index, const vector<double> &data,
                              bool entropy);

double loglikelihood(int start, int end, const vector<vector<int>> &index, bool entropy);


template<typename T>
void print_out_table(const vector<vector<T>>& vec, int width = 8, int precision = 3);
#endif //GOCTEST_CUTTING_CLUSTER_DP_COMPRESSED_H
