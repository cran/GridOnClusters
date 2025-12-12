//
// Created by Jiandong Wang on 25-2-6.
//

#include "utils.h"

#include <iostream>

using namespace std;

double getMean(vector<double> data)
{
    const size_t size = data.size();
    if (size == 0)
    {
        throw runtime_error("Cannot compute mean of an empty vector.");
    }

    const double sum = accumulate(data.begin(), data.end(), 0.0);
    const double mean = sum / size;

    return mean;
}

// Function to compute the column-wise mean for high-dimensional data
vector<double> computeMeans(const vector<vector<double>>& data)
{
    size_t numPoints = data.size();
    if (numPoints == 0) return {}; // Handle empty case

    size_t dim = data[0].size();
    vector<double> means(dim);

    for (size_t j = 0; j < dim; ++j)
    {
        vector<double> column(numPoints);
        for (size_t i = 0; i < numPoints; ++i)
        {
            column[i] = data[i][j]; // Extract j-th dimension values
        }
        means[j] = getMean(column);
    }
    return means;
}

double getMedian(vector<double> data)
{
    const size_t size = data.size();
    if (size == 0)
    {
        throw runtime_error("Cannot compute median of an empty vector.");
    }

    nth_element(data.begin(), data.begin() + size / 2, data.end());
    double median = data[size / 2];

    if (size % 2 == 0)
    {
        nth_element(data.begin(), data.begin() + size / 2 - 1, data.end());
        median = (median + data[size / 2 - 1]) / 2.0;
    }

    return median;
}

// Function to compute the column-wise median for high-dimensional data
vector<double> computeMedians(const vector<vector<double>>& data)
{
    size_t numPoints = data.size();
    if (numPoints == 0) return {}; // Handle empty case

    size_t dim = data[0].size();
    vector<double> medians(dim);

    for (size_t j = 0; j < dim; ++j)
    {
        vector<double> column(numPoints);
        for (size_t i = 0; i < numPoints; ++i)
        {
            column[i] = data[i][j]; // Extract j-th dimension values
        }
        medians[j] = getMedian(column);
    }
    return medians;
}

// Function to compute the column-wise centrals for high-dimensional data
vector<double> computeCentrals(const vector<vector<double>>& data, const string& method)
{
    size_t numPoints = data.size();
    if (numPoints == 0) return {}; // Handle empty case

    size_t dim = data[0].size();
    vector<double> centrals(dim);

    for (size_t j = 0; j < dim; ++j)
    {
        vector<double> column(numPoints);
        for (size_t i = 0; i < numPoints; ++i)
        {
            column[i] = data[i][j]; // Extract j-th dimension values
        }

        if (method == "median") centrals[j] = getMedian(column);
        else if (method == "mean") centrals[j] = getMean(column);
    }
    return centrals;
}

// Function to compute Euclidean distance to the centrals
vector<double> computeDistances(const vector<vector<double>>& data, const vector<double>& centrals)
{
    vector<double> distances(data.size());

    for (size_t i = 0; i < data.size(); ++i)
    {
        double sum = 0.0;
        for (size_t j = 0; j < centrals.size(); ++j)
        {
            double diff = data[i][j] - centrals[j];
            sum += diff * diff; // Squared difference
        }
        distances[i] = sqrt(sum); // Euclidean distance
    }
    return distances;
}
