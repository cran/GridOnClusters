//
// Created by Jiandong Wang on 25-2-6.
//

#ifndef UTILS_H
#define UTILS_H
#include <vector>
#include <string>
#include <algorithm>
#include <numeric> // For std::accumulate
#include <cmath>
#include <iostream>
#include <stdexcept>

template <typename T>
std::vector<std::vector<T>> getSubsetByIndices(
    const std::vector<std::vector<T>>& data,
    const std::vector<int>& indices)
{
    std::vector<std::vector<T>> subset;
    subset.reserve(indices.size());

    for (int index : indices) {
        if (index >= 0 && index < static_cast<int>(data.size())) {
            subset.push_back(data[index]);
        }
    }
    return subset;
}

double getMean(std::vector<double> data);

std::vector<double> computeMeans(const std::vector<std::vector<double>>& data);

double getMedian(std::vector<double> data);

std::vector<double> computeMedians(const std::vector<std::vector<double>>& data);

std::vector<double> computeCentrals(const std::vector<std::vector<double>>& data, const std::string& method);

std::vector<double> computeDistances(const std::vector<std::vector<double>>& data, const std::vector<double>& centrals);

#endif //UTILS_H
