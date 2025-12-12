//
// Created by Jiandong Wang on 2/15/20.
//
// Copyright (c) NMSU Song lab


#ifndef JOINT_GRID_CLUSTERS_H
#define JOINT_GRID_CLUSTERS_H


#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <numeric>
#include <cmath>


using namespace std;

const int DEBUG = 1;

struct grid {
    vector<vector<double> > lines;
    vector<vector<double>> medians;
    vector<vector<double> > median_distance;
    vector<vector<double> > mean_distance;
    vector<double> distance;
    vector<double> dimensional_bic;
    vector<size_t> num_lines;
    vector<vector<size_t>> conti_table;
    vector<vector<double> > lines_after_removed;
    vector<size_t> num_lines_after;
    vector<bool> line_removed;
    vector<vector<size_t>> conti_table_after;
    vector<vector<int>> discr_data;
    vector<vector<int>> discr_data_after;
    double original_BIC;
    double removed_BIC;
    double original_purity;
    double removed_purity;
    double upsilon;
};

struct Point {
    double point;
    int label;
};

// static const vector<double> *spx_c;

class Cluster {
private:
    unsigned int num_clusters = 0;
    unsigned int sample_size = 0;
    unsigned int sample_dims = 0;
    vector<vector<vector<double> > > cluster_points;
    vector<vector<double>> order;
    // cluster_label< dimension< points > >
    //vector<vector<double> > cluster_centers;
    vector<vector<double> > cluster_medians;
    vector<vector<double>> points;
    //  dimension< points >
    vector<int> labels;
    bool has_points = false;
    bool split = false;
    bool has_medians = false;
    bool has_ordered = false;
    grid grids;

    Cluster(vector<vector<double> > medians);

    Cluster(vector<vector<double> > medians, vector<vector<vector<double> > > data);

public:
    Cluster(vector<int> labels, vector<vector<double> > data);

    Cluster(vector<int> labels, vector<vector<double> > medians, vector<vector<double> > data);

    Cluster(int k, vector<int> labels, vector<vector<double> > data);

    void split_data();

    void cal_medians();

    vector<vector<double> > get_all_points();

    vector<double> get_all_points_by_dim(int dim);

    vector<vector<double>> convert_all_points_by_dim_cluster();

    vector<vector<double>> cal_order();

    vector<vector<vector<double> > > get_cluster_points();

    vector<vector<double> > get_points(int dim);

    vector<double> get_points(int label, int dim);

    double get_points(int label, int dim, int index);

    vector<vector<double> > get_medians();

    vector<double> get_medians(int label);

    vector<vector<double>> get_order();

    vector<int> get_labels();

    int get_K();

    void set_grids(grid G);

    grid get_grids();

    vector<int> sort_clusters(int dim);

    int get_dims();

    int get_sample_size();

    bool has_all_points();

    bool split_by_cluster();

    bool has_median();
};


#endif //JOINT_GRID_CLUSTERS_H
