// Joint_Grid.cpp
//
// Created by Jiandong Wang on 2/15/20.
//
// Copyright (c) NMSU Song lab


#include "Joint_Grid.h"

/*
vector<vector<double> > Find_Grid(Cluster &clusters) {
    int dims = clusters.get_dims();
    vector<vector<double >> lines_multi_dim;
    vector<int> num_lines;
    for (int i = 1; i < dims + 1; ++i) {
        vector<double> lines = Find_1D_Grid(clusters, i);
        lines_multi_dim.push_back(lines);
        num_lines.push_back(lines.size());
    }
    grid G;
    G.lines = lines_multi_dim;
    G.num_lines = num_lines;
    clusters.set_grids(G);
    return lines_multi_dim;
}
 */

vector<vector<double> > Find_Grid(Cluster & clusters) {
    int dims = clusters.get_dims();
    vector<vector<double >> lines_multi_dim(dims);
    vector<int> num_lines(dims);
    for (int i = 0; i < dims; ++i) {
        vector<double> lines = Find_1D_Grid(clusters, i+1);
        lines_multi_dim[i] = lines;
        num_lines[i] = lines.size();
    }
    grid G;
    G.lines = lines_multi_dim;
    G.num_lines = num_lines;
    clusters.set_grids(G);
    return lines_multi_dim;
}

vector<double> Find_1D_Grid(Cluster &clusters, int dim_input) {
    int dim = dim_input - 1;
    vector<int> order = clusters.sort_clusters(dim);
    vector<double> c1, c2;
    double center_1, center_2;
    vector<double> lines;
    vector<double> overlap_lines;
    vector<double> err_sum;
    for (int i = 0; i < order.size() - 1; ++i) {
        // extract two clusters on that dimension
        c1 = clusters.get_points(order[i])[dim];
        c2 = clusters.get_points(order[i + 1])[dim];
        center_1 = clusters.get_centers(order[i])[dim];
        center_2 = clusters.get_centers(order[i + 1])[dim];

        vector<vector<double> > c_index;

        c_index = prep_index(c1, c2, center_1, center_2);


        //searching
        int len_c1 = c1.size();
        int len_c2 = c2.size();
        bool overlap = false;

        double line;
        if (c_index[0].size() == 0) {
            overlap = true;
            line = (center_1 + center_2) / 2.0;
        } else {
            line = binary_search_index(c_index, 0, c_index[0].size() - 1, len_c1, len_c2, overlap, err_sum);
        }
        if (not overlap) {
            lines.push_back(line);
        } else if (overlap) {
            overlap_lines.push_back(line);
        }
    }
    if (lines.size() == 0) {
        auto posi = min_element(err_sum.begin(), err_sum.end());
        lines.push_back(overlap_lines[posi - err_sum.begin()]);
    }

    return lines;

}

vector<vector<double> > prep_index(vector<double> &c1, vector<double> &c2, double center_1, double center_2) {
    // sort two clusters
    sort(c1.begin(), c1.end());
    sort(c2.begin(), c2.end());
    auto c1_g = lower_bound(c1.begin(), c1.end(), center_1);
    auto c1_l = upper_bound(c1.begin(), c1.end(), center_2);
    vector<double> c1_new(c1_g, c1_l);

    auto c2_g = lower_bound(c2.begin(), c2.end(), center_1);
    auto c2_l = upper_bound(c2.begin(), c2.end(), center_2);
    vector<double> c2_new(c2_g, c2_l);


    // int len_c1 = c1_new.size();  Commented MS. March 17, 2020
    // int len_c2 = c2_new.size();  Commented MS. March 17, 2020

    // merge them together
    vector<double>::iterator iter_c1 = c1_new.begin();
    vector<double>::iterator iter_c2 = c2_new.begin();
    vector<vector<double> > c_all = vector<vector<double> >(3, vector<double>(c1_new.size() + c2_new.size(), 0));

    vector<double>::iterator iter_c_all_0 = c_all[0].begin();
    vector<double>::iterator iter_c_all_1 = c_all[1].begin();
    vector<double>::iterator iter_c_all_2 = c_all[2].begin();
    int prev_1 = c1.end() - c1_g + 1;
    int prev_2 = 0;
    for (; iter_c_all_0 != c_all[0].end(); *iter_c_all_0++, *iter_c_all_1++, *iter_c_all_2++) {

        if (iter_c1 == c1_new.end()) {//c1 is empty, put c2
            *iter_c_all_0 = *iter_c2;
            iter_c2++;
            prev_2++;
            prev_1 = c1.end() - c1_l;

        } else if (iter_c2 == c2_new.end()) {//c2 is empty, put c1
            *iter_c_all_0 = *iter_c1;
            prev_1--;
            iter_c1++;
            prev_2 = c2_g - c2.begin();

        } else if (*iter_c1 <= *iter_c2) {//put c1
            *iter_c_all_0 = *iter_c1;
            prev_1--;
            iter_c1++;

        } else if (*iter_c1 > *iter_c2) {//put c2
            *iter_c_all_0 = *iter_c2;
            iter_c2++;
            prev_2++;
        } else {
            throw "Condition Error!!!\n";
        }

        *iter_c_all_1 = prev_1;
        *iter_c_all_2 = prev_2;
    }
    return c_all;
}

double binary_search_index(const vector<vector<double> > &c_index, const int left, const int right, const int size_c1,
                           const int size_c2, bool &overlap, vector<double> &err_sum) {
    int line_index = 0;
    int line_index_left = 0;
    int line_index_right = 0;
    double line = 0;
    int index_c1 = 0;
    int index_c2 = 0;
    bool even = false;

    if (right < left) {
        throw "VALUE ERROR!\nRight bound less than left bound!\n";
    }
    if ((right - left) % 2 == 0) {
        line_index = left + (right - left) / 2;
        index_c1 = c_index[1][line_index];
        index_c2 = c_index[2][line_index];
        line = c_index[0][line_index];
    } else if ((right - left) % 2 == 1) {
        line_index_left = left + (right - left) / 2;
        line_index_right = left + (right - left) / 2 + 1;
        index_c1 = c_index[1][line_index_right];
        index_c2 = c_index[2][line_index_left];
        line = (c_index[0][line_index_left] + c_index[0][line_index_right]) / 2.0;
        even = true;
    }

    double e_rate_c1 = index_c1 / float(size_c1);
    double e_rate_c2 = index_c2 / float(size_c2);

    if (abs(e_rate_c1 - e_rate_c2) <= 0.005 and e_rate_c1 + e_rate_c2 >= 0.50) {
        overlap = true;
        err_sum.push_back(e_rate_c1 + e_rate_c2);
    } else if (abs(e_rate_c1 - e_rate_c2) <= 0.005 and e_rate_c1 + e_rate_c2 < 0.50) {
        overlap = false;
    } else if (left == right) {
        if (e_rate_c1 + e_rate_c2 >= 0.50) {
            overlap = true;
            err_sum.push_back(e_rate_c1 + e_rate_c2);
        } else { overlap = false; }
    } else if (e_rate_c1 > e_rate_c2) {
        // moving right
        if (even) {
            line_index = line_index_right;
        }
        return binary_search_index(c_index, line_index, right, size_c1, size_c2, overlap, err_sum);
    } else if (e_rate_c1 < e_rate_c2) {
        //moving left
        if (even) {
            line_index = line_index_left;
        }
        return binary_search_index(c_index, left, line_index, size_c1, size_c2, overlap, err_sum);
    }
    return line;
}






















