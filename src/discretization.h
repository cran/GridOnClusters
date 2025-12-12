//
// Created by Jiandong Wang on 2022/11/12.
//

#ifndef GOCTEST_DISCRETIZATION_H
#define GOCTEST_DISCRETIZATION_H

#include <vector>
#include <thread>
#include <mutex>
#include "Clusters.h"

using namespace std;

class my_compare {
public:
    static vector<vector<double>> data;
    vector<size_t> order;
    //static const vector<double> *spx;
    static size_t dim;

    my_compare(vector<vector<double>> input_data, vector<size_t> input_order) {
        this->data = input_data;
        this->order = input_order;
    }

    static bool compi(size_t i, size_t j) {
        return data[dim][i] < data[dim][j];
    }

    void sort_on_dim(size_t dim) {
        this->dim = dim;
        sort(this->order.begin(), this->order.end(), this->compi);
    }
};

vector<vector<int>>
discretize_data(Cluster &cluster, vector<vector<double>> lines_all_dim, size_t n_thread, size_t step,
                bool skip = false);

void discretization_data_thread(vector<vector<int>> &discr_data, const vector<vector<double>> &data,
                                const vector<vector<double>> &order, const vector<vector<double>> &lines_all_dim,
                                size_t nsample, size_t dim_u, size_t dim_l, bool skip);

void _test_discretize_data();

#endif //GOCTEST_DISCRETIZATION_H
