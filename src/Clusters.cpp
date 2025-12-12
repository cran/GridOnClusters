//
// Created by Jiandong Wang on 2/15/20.
//
// Copyright (c) NMSU Song lab

#include "Clusters.h"

static const vector<double> *spx_c;

// testing purpose only
Cluster::Cluster(vector<vector<double> > medians) {
    this->num_clusters = medians.size();
    this->cluster_medians = medians;
}


// testing purpose only
Cluster::Cluster(vector<vector<double> > medians, vector<vector<vector<double> > > data) {
    this->num_clusters = medians.size();
    this->cluster_medians = medians;
    this->cluster_points = data;
}

Cluster::Cluster(vector<int> labels, vector<vector<double> > data) {
    this->num_clusters = *max_element(labels.begin(), labels.end());
    this->sample_size = data.size();
    this->sample_dims = data[0].size();
    this->labels = labels;
    this->points = data;
    this->has_points = true;
    this->split = false;
    this->has_medians = false;
}

Cluster::Cluster(int k, vector<int> labels, vector<vector<double> > data) {
    this->num_clusters = k;
    this->sample_size = data.size();
    this->sample_dims = data[0].size();
    this->labels = labels;
    this->points = data;
    this->has_points = true;
    this->split = false;
    this->has_medians = false;
}

Cluster::Cluster(vector<int> labels, vector<vector<double> > medians, vector<vector<double> > data) {
    this->points = data;
    this->cluster_medians = medians;
    this->labels = labels;
    this->num_clusters = medians.size();
    this->sample_dims = data[0].size();
    this->has_points = true;
    this->has_medians = true;
    this->split = false;
}

void Cluster::split_data() {
    vector<vector<vector<double>>> c_data(num_clusters,
                                          vector<vector<double>>(sample_dims, vector<double>(sample_size, 0)));
    vector<size_t> clusters_size(num_clusters, 0);
    vector<int>::iterator label_iter = this->labels.begin();
    vector<vector<double> >::iterator data_iter = this->points.begin();
    int label;
    //loop for clusters
    for (; label_iter != this->labels.end() and data_iter != this->points.end(); ++label_iter, ++data_iter) {
        //loop for dimensions
        label = *label_iter;
        for (size_t i = 0; i < this->sample_dims; ++i) {
            c_data[label][i][clusters_size[label]] = (*data_iter)[i];
        }
        clusters_size[label]++;
    }

    for (size_t i = 0; i < this->num_clusters; ++i) {
        for (size_t j = 0; j < this->sample_dims; ++j) {
            (c_data[i][j]).resize(clusters_size[i]);
        }
    }

    this->cluster_points = c_data;
}

void Cluster::cal_medians() {
    this->cluster_medians = vector<vector<double> >(this->num_clusters, vector<double>(sample_dims, 0));
    double mid;
    int mid_index;
    for (size_t i = 0; i < this->num_clusters; ++i) {
        for (size_t j = 0; j < sample_dims; ++j) {
            vector<double> data = this->cluster_points[i][j];
            sort(data.begin(), data.end());
            if (data.size() % 2 == 0) {//mid_point is in between
                mid_index = data.size() / 2;
                mid = (data[mid_index] + data[mid_index - 1]) / 2.0;
            } else {//mid_point is in data
                mid_index = data.size() / 2.0;
                mid = data[floor(mid_index)];
            }
            this->cluster_medians[i][j] = mid;
        }
    }
}

vector<vector<vector<double> > > Cluster::get_cluster_points() {
    return this->cluster_points;
}

vector<vector<double> > Cluster::get_all_points() {
    return this->points;
}

vector<double> Cluster::get_all_points_by_dim(int dim) {
    vector<double> points_dim(this->sample_size);
    auto dim_iter = points_dim.begin();
    for (auto iter = this->points.begin(); iter != this->points.end(); ++iter, ++dim_iter) {
        (*dim_iter) = (*iter)[dim];
    }
    return points_dim;
}

vector<vector<double>> Cluster::convert_all_points_by_dim_cluster() {
    vector<vector<double>> points_by_dim(this->sample_dims, vector<double>(sample_size));
    vector<int> index(this->sample_dims, 0);

    for (size_t i = 0; i < this->sample_size; ++i) {
        for (size_t j = 0; j < this->sample_dims; ++j) {
            points_by_dim[j][i] = this->points[i][j];
        }
    }

    return points_by_dim;
}

static bool compi_c(size_t i, size_t j) {
    return (*spx_c)[i] < (*spx_c)[j];
}

vector<vector<double>> Cluster::cal_order() {
    auto data_by_dim = this->convert_all_points_by_dim_cluster();
    if (!this->has_ordered) {
        this->order = vector<vector<double>>(sample_dims, vector<double>(sample_size, 0));
        vector<double> order_sample(this->sample_size);
        iota(order_sample.begin(), order_sample.end(), 0);
        for (size_t dim = 0; dim < sample_dims; ++dim) {
            spx_c = &data_by_dim[dim];
            auto data_order = order_sample;
            sort(data_order.begin(), data_order.end(), compi_c);
            this->order[dim] = std::move(data_order);
        }
        this->has_ordered = true;
    }
    return data_by_dim;
}

vector<vector<double> > Cluster::get_points(int dim) {
    vector<vector<double>> points_dim(this->num_clusters);
    auto iter_res = points_dim.begin();
    for (auto iter = this->cluster_points.begin(); iter != this->cluster_points.end(); ++iter, ++iter_res) {
        (*iter_res) = (*iter)[dim];
    }
    return points_dim;
}

vector<double> Cluster::get_points(int label, int dim) {
    return this->cluster_points[label][dim];
}

double Cluster::get_points(int label, int dim, int index) {
    return this->cluster_points[label][dim][index];
}

vector<vector<double> > Cluster::get_medians() {
    return this->cluster_medians;
}

vector<double> Cluster::get_medians(int label) {
    return this->cluster_medians[label];
}

vector<vector<double>> Cluster::get_order() {
    return this->order;
}

int Cluster::get_K() {
    return this->num_clusters;
};

void Cluster::set_grids(grid G) {
    this->grids = G;
}

vector<int> Cluster::get_labels() {
    return this->labels;
}

grid Cluster::get_grids() {
    return this->grids;
}

vector<int> Cluster::sort_clusters(int dim) {
    /// sort clusters by their median
    // initialize original index locations
    vector<int> idx(this->num_clusters);
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [this, dim](size_t i1, size_t i2) { return this->cluster_medians[i1][dim] < this->cluster_medians[i2][dim]; });
    return idx;
}

int Cluster::get_dims() {
    return this->sample_dims;
}

int Cluster::get_sample_size() { return this->sample_size; }

bool Cluster::has_all_points() {
    return this->has_points;
}

bool Cluster::split_by_cluster() {
    return this->split;
}

bool Cluster::has_median() {
    return this->has_medians;
}
