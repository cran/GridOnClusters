//
// Created by Jiandong Wang on 2025/2/4.
//

#include <R.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "evaluation.h"
#include "MyHash.h"
#include "utils.h"

using namespace std;

vector<vector<double>>  evaluation_Distances(const vector<vector<double>>& data, const MyHash& hct, const string& method)
{
    vector<vector<double>>  distances=vector<vector<double>>(hct.m_table.size());
    auto iter_dist = distances.begin();
    for (auto p = hct.m_table.begin(); p != hct.m_table.end(); ++p)
    {
        auto data_sub = getSubsetByIndices(data, p->second.members);
        vector<double> centrals;
        if (method == "median") centrals = computeMedians(data_sub);
        else if (method == "mean") centrals = computeMeans(data_sub);
        auto temp_dist = computeDistances(data_sub, centrals);
        *iter_dist = temp_dist;
        ++iter_dist;
    }
    return distances;
}

void test_mid_dist()
{
    // Example high-dimensional dataset (5 points, 3 dimensions)
    vector<vector<double>> data = {
        {2.0, 3.0, 4.0},
        {5.0, 6.0, 7.0},
        {1.0, 1.5, 2.0},
        {3.0, 3.5, 4.5},
        {4.0, 4.5, 5.0}
    };

    // Compute median for each dimension
    vector<double> median = computeMedians(data);

    // Compute distances from median
    vector<double> distances = computeDistances(data, median);

    // Output results
    Rprintf("Median of each dimension: ");
    for (double m : median) {
      // cout << m << " "; 
      Rprintf("%f ", m);
    }
    
    Rprintf("\nEuclidean distances: ");
    for (double d : distances) {
      // cout << d << " "; 
      Rprintf("%f ", d);
    }
    Rprintf("\n");
}
