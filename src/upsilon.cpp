//
// Created by Jiandong Wang on 2022/2/2.
//

#include "upsilon.h"

// double upsilon_stat(const vector<vector<unsigned long> > &data) {
//     double upsilon = 0.0;
//     if (data.empty()) {
//         return upsilon;
//     } else if (data.empty()) {
//         return upsilon;
//     }
// 
//     int n = 0;
//     std::vector<int> colsums((int) data[0].size(), 0);
//     std::vector<int> rowsums((int) data.size(), 0);
// 
//     for (size_t i = 0; i < data.size(); i++) {
//         for (size_t j = 0; j < data[i].size(); j++) {
//             n += data[i][j];
//             colsums[j] += data[i][j];
//             rowsums[i] += data[i][j];
//         }
//     }
// 
//     if (n == 0)return upsilon;
// 
//     size_t nrows = data.size();  // number of rows
//     size_t ncols = data[0].size();  // number of columns
// 
//     for (size_t i = 0; i < nrows; ++i) {
//         // Expected cound for cell (i,j):
//         for (size_t j = 0; j < ncols; ++j) {
//             double eij = rowsums[i] * colsums[j] / (double) n;
//             upsilon += (data[i][j] - eij) * (data[i][j] - eij) * nrows * ncols / n;
// 
//         }
//     }
//     return upsilon;
// }
// 
// double upsilon_stat(const vector<vector<int>> &data) {
//     double upsilon = 0.0;
//     if (data.empty()) {
//         return upsilon;
//     } else if (data.empty()) {
//         return upsilon;
//     }
// 
//     int n = 0;
//     std::vector<int> colsums((int) data[0].size(), 0);
//     std::vector<int> rowsums((int) data.size(), 0);
// 
//     for (size_t i = 0; i < data.size(); i++) {
//         for (size_t j = 0; j < data[i].size(); j++) {
//             n += data[i][j];
//             colsums[j] += data[i][j];
//             rowsums[i] += data[i][j];
//         }
//     }
// 
//     if (n == 0)return upsilon;
// 
//     size_t nrows = data.size();  // number of rows
//     size_t ncols = data[0].size();  // number of columns
// 
//     for (size_t i = 0; i < nrows; ++i) {
//         // Expected cound for cell (i,j):
//         for (size_t j = 0; j < ncols; ++j) {
//             double eij = rowsums[i] * colsums[j] / (double) n;
//             upsilon += (data[i][j] - eij) * (data[i][j] - eij) * nrows * ncols / n;
// 
//         }
//     }
//     return upsilon;
// };

double upsilon_stat(const vector<vector<unsigned long>> &data, double &df, double effect) {
    vector<vector<int>> data_int(data.size(), vector<int>(data[0].size(), 0));
    //copy(data.begin(), data.end(), data_int);

    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[0].size(); ++j) {
            data_int[i][j] = (int) data[i][j];
        }
    }

    double stat = 0;
    double ef;
    upsilon_test_by_hashing(data_int, stat, df, ef);
    return stat;
};

double upsilon_stat(const vector<vector<int>> &data, double &df, double effect) {
    double ef;
    double stat = 0;
    upsilon_test_by_hashing(data, stat, df, ef);
    return stat;
};
