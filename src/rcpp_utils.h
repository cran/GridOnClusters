//
// Created by Jiandong Wang on 25-2-6.
//

#ifndef RCPP_UTILS_H
#define RCPP_UTILS_H

#include <vector>
#include <Rcpp.h>

// Function to convert a 2D vector into a Rcpp::List
template <typename T>
Rcpp::List convertToRList(const std::vector<std::vector<T>>& table) {
    Rcpp::List table_o(table.size()); // Preallocate memory

    // Use std::transform for efficient mapping
    std::transform(table.begin(), table.end(), table_o.begin(),
                   [](const std::vector<T>& row) { return Rcpp::wrap(row); });

    return table_o;
}

template <typename T>
std::vector<std::vector<T> > convertToVector(const Rcpp::NumericMatrix & data) {
    std::vector<std::vector<T> > c_data(data.nrow(), std::vector<T>(data.ncol(), 0));

    for (auto i = 0; i < data.nrow(); ++i) {
        for (auto j = 0; j < data.ncol(); ++j) {
            c_data[i][j] = data(i, j);
        }
    }
    return c_data;
}

#endif //RCPP_UTILS_H
