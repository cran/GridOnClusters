//
// Created by Jiandong Wang on 2022/12/7.
//

#ifndef GOCTEST_EVALUATION_H
#define GOCTEST_EVALUATION_H

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <algorithm>
#include "MyHash.h"

double evaluation_entropy(double N, double Q, double K, std::vector<std::vector<double>> S,
                          std::vector<std::vector<double>> S_entropy);

double evaluation_likelihood(double N, double Q, double K, std::vector<std::vector<double>> S);

double evaluation_counts(const std::vector<double>& points, const std::vector<int>& labels,
                         std::vector<std::vector<int>> index, std::vector<int> posi,
                         std::vector<std::vector<int>>& table);

void shifted_data_variance(const std::vector<double>& x, int left, int right, double& mean, double& variance);

std::vector<std::vector<double>> split_into_channels(const std::vector<double>& points, const std::vector<int>& labels,
                                                     unsigned int K, unsigned int N);

double purity(std::vector<std::vector<size_t>> table, size_t N);

std::vector<std::vector<double>> evaluation_Distances(const std::vector<std::vector<double>>& data, const MyHash& hct, const std::string& method);
#endif //GOCTEST_EVALUATION_H
