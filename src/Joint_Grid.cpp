//
// Created by Jiandong Wang on 2/15/20.
//
// Copyright (c) NMSU Song lab


#include "Joint_Grid.h"

vector<vector<double> > Find_Grid(Cluster &clusters) {
    int dims = clusters.get_dims();
    vector<vector<double> > lines_multi_dim(dims);
    vector<int> num_lines(dims);
    for (int i = 0; i < dims; ++i) {
        vector<double> lines = Find_1D_Grid(clusters, i + 1);
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
    double mid_1, mid_2;
    vector<double> lines(order.size(), ULONG_MAX);
    vector<double> overlap_lines;
    vector<double> err_sum;
    bool line_change = false;
    for (size_t i = 0; i < int(order.size()) - 1; ++i) {
        // extract two clusters on that dimension
        c1 = clusters.get_points(order[i])[dim];
        c2 = clusters.get_points(order[i + 1])[dim];
        mid_1 = clusters.get_medians(order[i])[dim];
        mid_2 = clusters.get_medians(order[i + 1])[dim];

        vector<vector<double> > c_index;

        c_index = prep_index(c1, c2, mid_1, mid_2);


        //searching
        int len_c1 = c1.size();
        int len_c2 = c2.size();
        bool overlap = false;

        double line;
        if (c_index[0].size() == 0) {
            overlap = true;
            line = (mid_1 + mid_2) / 2.0;
        } else {
            line = binary_search_index(c_index, 0, c_index[0].size() - 1, len_c1, len_c2, overlap, err_sum);
        }

        if (not overlap) {
            lines[i]=(line);
            line_change = true;
        } else if (overlap) {
            overlap_lines.push_back(line);
        }
    }
    if (not line_change) {
        // for (size_t n = 0; n < err_sum.size(); ++n) {
        //     cout << "line: " << overlap_lines[n] << " error sum: " << err_sum[n] << endl;
        // }
        auto posi = min_element(err_sum.begin(), err_sum.end());
        lines[0]=(overlap_lines[posi - err_sum.begin()]);
        //cout << "Final line is: " << lines[0] << endl;
    }
    return lines;

}

vector<vector<double> > prep_index(vector<double> &c1, vector<double> &c2, double median_1, double median_2) {
    // sort two clusters
    sort(c1.begin(), c1.end());
    sort(c2.begin(), c2.end());
    auto c1_l = lower_bound(c1.begin(), c1.end(), median_1);
    auto c1_r = upper_bound(c1.begin(), c1.end(), median_2);
    vector<double> c1_new(c1_l, c1_r);

    auto c2_l = lower_bound(c2.begin(), c2.end(), median_1);
    auto c2_r = upper_bound(c2.begin(), c2.end(), median_2);
    vector<double> c2_new(c2_l, c2_r);

    // merge them together
    vector<double>::iterator iter_c1 = c1_new.begin();
    vector<double>::iterator iter_c2 = c2_new.begin();
    vector<vector<double> > data_index = vector<vector<double> >(3, vector<double>(c1_new.size() + c2_new.size(), 0));

    vector<double>::iterator real_data = data_index[0].begin();
    vector<double>::iterator c1_index = data_index[1].begin();
    vector<double>::iterator c2_index = data_index[2].begin();
    int prev_1 = distance(c1_l, c1.end());
    int prev_2 = distance(c2.begin(), c2_l);
    for (; real_data != data_index[0].end(); *real_data++, *c1_index++, *c2_index++) {

        if (iter_c1 == c1_new.end()) {//c1 is empty, put c2
            *real_data = *iter_c2;
            iter_c2++;
            prev_2++;
            prev_1 = distance(c1_r, c1.end());

        } else if (iter_c2 == c2_new.end()) {//c2 is empty, put c1
            *real_data = *iter_c1;
            prev_1--;
            iter_c1++;
            prev_2 = distance(c2.begin(), c2_r - 1);

        } else if (*iter_c1 <= *iter_c2) {//put c1
            *real_data = *iter_c1;
            prev_1--;
            iter_c1++;

        } else if (*iter_c1 > *iter_c2) {//put c2
            *real_data = *iter_c2;
            iter_c2++;
            prev_2++;
        } else {
            throw "Condition Error!!!\n";
        }

        *c1_index = prev_1;
        *c2_index = prev_2;
    }

    if (data_index[0].size() == 0) {
        vector<vector<double> > mean_index = vector<vector<double> >(3, vector<double>(0, 0));
        return mean_index;
    }
    vector<vector<double> > mean_index = vector<vector<double> >(3, vector<double>(data_index[0].size() - 1, 0));

    for (int i = 0; i < int(data_index[0].size()) - 1; ++i) {
        double mean = (data_index[0][i] + data_index[0][i + 1]) / 2.0;
        int c1_index = data_index[1][i];
        int c2_index = data_index[2][i];
        mean_index[0][i] = mean;
        mean_index[1][i] = c1_index;
        mean_index[2][i] = c2_index;
    }

    return mean_index;
}

double binary_search_index(const vector<vector<double> > &c_index, const int left, const int right, const int size_c1,
                           const int size_c2, bool &overlap, vector<double> &err_sum) {
//    int line_index = 0;
//    int line_index_left = 0;
//    int line_index_right = 0;
//    double line = 0;
//    int index_c1 = 0;
//    int index_c2 = 0;
//    bool even = false;

    if (right < left) {
        throw "VALUE ERROR!\nRight bound less than left bound!\n";
    }
//    if ((right - left) % 2 == 0) {
//        line_index = left + (right - left) / 2;
//        index_c1 = c_index[1][line_index];
//        index_c2 = c_index[2][line_index];
//        line = c_index[0][line_index];
//    } else if ((right - left) % 2 == 1) {
//        line_index_left = left + (right - left) / 2;
//        line_index_right = left + (right - left) / 2 + 1;
//        index_c1 = c_index[1][line_index_right];
//        index_c2 = c_index[2][line_index_left];
//        line = (c_index[0][line_index_left] + c_index[0][line_index_right]) / 2.0;
//        even = true;
//    }

    int line_index = ceil(left + (right - left) / 2.0);
    int c1_index = c_index[1][line_index];
    int c2_index = c_index[2][line_index];
    double e_rate_c1 = c1_index / float(size_c1);
    double e_rate_c2 = c2_index / float(size_c2);

    if (abs(e_rate_c1 - e_rate_c2) <= 0.005 and e_rate_c1 + e_rate_c2 >= 0.50) {
        overlap = true;
        err_sum.push_back(e_rate_c1 + e_rate_c2);
    } else if (abs(e_rate_c1 - e_rate_c2) <= 0.005 and e_rate_c1 + e_rate_c2 < 0.50) {
        overlap = false;
    } else if (left == right - 1) {
        double left_e_rate_c1 = c_index[1][line_index - 1] / float(size_c1);
        double left_e_rate_c2 = c_index[2][line_index - 1] / float(size_c2);

        if ((left_e_rate_c1 + left_e_rate_c2) < (e_rate_c1 + e_rate_c2)) {
            line_index = line_index - 1;
            e_rate_c1 = left_e_rate_c1;
            e_rate_c2 = left_e_rate_c2;
        }
        if (e_rate_c1 + e_rate_c2 >= 0.50) {
            overlap = true;
            err_sum.push_back(e_rate_c1 + e_rate_c2);
        } else { overlap = false; }
    } else if (left == right) {
        if (e_rate_c1 + e_rate_c2 >= 0.50) {
            overlap = true;
            err_sum.push_back(e_rate_c1 + e_rate_c2);
        } else { overlap = false; }
    } else if (e_rate_c1 > e_rate_c2) {
        // moving right
//        if (even) {
//            line_index = line_index_right;
//        }
        return binary_search_index(c_index, line_index, right, size_c1, size_c2, overlap, err_sum);
    } else if (e_rate_c1 < e_rate_c2) {
        //moving left
//        if (even) {
//            line_index = line_index_left;
//        }
        return binary_search_index(c_index, left, line_index, size_c1, size_c2, overlap, err_sum);
    }
    return c_index[0][line_index];
}






















