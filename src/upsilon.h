//
// Created by Jiandong Wang on 2022/2/2.
//

#ifndef GOCTEST_UPSILON_H
#define GOCTEST_UPSILON_H

#include <vector>
#include "HashedUpsilon.h"

using namespace std;

//double upsilon_stat(const vector<vector<unsigned long>> &data);

//double upsilon_stat(const vector<vector<int>> &data);

#include <vector>
#include <iostream>

template <typename T>
double upsilon_stat(const std::vector<std::vector<T>>& data) {
   double upsilon = 0.0;
   if (data.empty()) {
      return upsilon;
   }
   
   int n = 0;
   std::vector<int> colsums(data[0].size(), 0);
   std::vector<int> rowsums(data.size(), 0);
   
   for (size_t i = 0; i < data.size(); i++) {
      for (size_t j = 0; j < data[i].size(); j++) {
         n += data[i][j];
         colsums[j] += data[i][j];
         rowsums[i] += data[i][j];
      }
   }
   
   if (n == 0) return upsilon;
   
   size_t nrows = data.size();  // number of rows
   size_t ncols = data[0].size();  // number of columns
   
   for (size_t i = 0; i < nrows; ++i) {
      // Expected count for cell (i,j):
      for (size_t j = 0; j < ncols; ++j) {
         double eij = rowsums[i] * colsums[j] / (double) n;
         upsilon += (data[i][j] - eij) * (data[i][j] - eij) * nrows * ncols / n;
      }
   }
   return upsilon;
}
   
double upsilon_stat(const vector<vector<unsigned long>> &data, double &df, double effect);

double upsilon_stat(const vector<vector<int>> &data, double &df, double effect);

#endif //GOCTEST_UPSILON_H
