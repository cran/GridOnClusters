//
// Created by Jiandong Wang on 2/15/20.
//
// Copyright (c) NMSU Song lab

#include "Clusters.h"

// testing purpose only
Cluster::Cluster(vector<vector<double> > centers) {
    this->num_clusters = centers.size();
    this->cluster_centers = centers;
}


// testing purpose only
Cluster::Cluster(vector<vector<double> > centers, vector<vector<vector<double> > > data) {
    this->num_clusters = centers.size();
    this->cluster_centers = centers;
    this->cluster_points = data;
}

Cluster::Cluster(vector<int> labels, vector<vector<double> > centers, vector<vector<double> > data) {
    // get number of clusters
    this->num_clusters = centers.size();
    int dims = centers[0].size();
    this->cluster_points = vector<vector<vector<double> > >(num_clusters,
                                                            vector<vector<double>>(dims, vector<double>(0)));

    // get clusters
    vector<int>::iterator label_iter = labels.begin();
    vector<vector<double> >::iterator data_iter = data.begin();
    //loop for clusters
    for (; label_iter != labels.end() and data_iter != data.end(); ++label_iter, ++data_iter) {
        //loop for dimensions
        for (int i = 0; i < dims; ++i) {
            this->cluster_points[*label_iter][i].push_back((*data_iter)[i]);
        }
    }
    this->cluster_centers = centers;
}

void Cluster::set_clusters(vector<int> labels, vector<vector<double> > centers, vector<vector<double> > data) {
    // get number of clusters
    this->num_clusters = centers.size();
    int dims = centers[0].size();
    this->cluster_points = vector<vector<vector<double> > >(num_clusters,
                                                            vector<vector<double>>(dims, vector<double>(0)));

    // get clusters
    vector<int>::iterator label_iter = labels.begin();
    vector<vector<double> >::iterator data_iter = data.begin();
    //loop for clusters
    for (; label_iter != labels.end() and data_iter != data.end(); ++label_iter, ++data_iter) {
        //loop for dimensions
        for (int i = 0; i < dims; ++i) {
            this->cluster_points[*label_iter][i].push_back((*data_iter)[i]);
        }
    }
    this->cluster_centers = centers;
}

vector<vector<vector<double> > > Cluster::get_points() {
    return this->cluster_points;
}

vector<vector<double> > Cluster::get_centers() {
    return this->cluster_centers;
}

void Cluster::set_grids(grid G) {
    this->grids = G;
}

grid Cluster::get_grids() {
    return this->grids;
}

vector<int> Cluster::sort_clusters(int dim) {
    // initialize original index locations
    vector<int> idx(this->num_clusters);
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [this, dim](size_t i1, size_t i2) { return this->cluster_centers[i1][dim] < this->cluster_centers[i2][dim]; });

    // if (DEBUG == 1) {
    //     vector<int>::iterator ite = idx.begin();
    //     for (; ite != idx.end(); ite++) {
    //         cout << *ite;
    //     }
    //     cout << endl;
    // }
    return idx;
}

int Cluster::get_dims() {
    return this->cluster_centers[0].size();
}

vector<vector<double> > Cluster::get_points(int index) {
    return this->cluster_points[index];
}


vector<double> Cluster::get_centers(int index) {
    return this->cluster_centers[index];
}
