//
// Created by Jiandong Wang on 2022/2/23.
//

#include "Cut_by_density.h"

vector<double> Cut_by_density_1D(Cluster &clusters, int dim_input, int min_bin_limit) {
    int dim = dim_input;// the current dimension, starting from 0
    vector<int> order = clusters.sort_clusters(dim);
    vector<double> lines(order.size(), 0);

    if (order.empty()) {
        return lines;
    }

    vector<double> c1, c2;
    double mid_1, mid_2;
    vector<double> overlap_lines;
    vector<double> err_sum;
    auto iter_lines = lines.begin();
    int num_bins = 1;
    int num_overlaps = 0;
    for (size_t i = 0; i + 1 < order.size(); ++i) {
        // extract two clusters on that dimension
        c1 = clusters.get_points(order[i], dim);
        c2 = clusters.get_points(order[i + 1], dim);
        mid_1 = clusters.get_medians(order[i])[dim];
        mid_2 = clusters.get_medians(order[i + 1])[dim];

        //check simple cluster
        bool c1_simple = false;
        bool c2_simple = false;
        // c1 is simple
        if (mid_1 == *min_element(c1.begin(), c1.end())) {
            c1_simple = true;
        }
        // c2 is simple
        if (mid_2 == *min_element(c2.begin(), c2.end())) {
            c2_simple = true;
        }
        // c1 and c2 both simple
        if (c1_simple and c2_simple) {

            // c1 == c2, go for overlap, set line by mid, skip rest
            if (mid_1 == mid_2) {
                overlap_lines.push_back(mid_1);
                num_overlaps++;
                continue;
            }
        }

        vector<vector<double> > c_index;

        c_index = prep_index(c1, c2, mid_1, mid_2);


        //searching
        int len_c1 = c1.size();
        int len_c2 = c2.size();
        bool overlap = false;

        double line;
        if (c_index[0].empty()) {
            overlap = true;
            line = (mid_1 + mid_2) / 2.0;
            err_sum.push_back(line);
        } else {
            line = binary_search_index(c_index, 0, c_index[0].size() - 1, len_c1, len_c2, overlap, err_sum);
        }

        if (not overlap) {
            *iter_lines = line;
            iter_lines++;
            num_bins++;
        } else if (overlap) {
            overlap_lines.push_back(line);
            num_overlaps++;
        }
    }
    if (num_bins < min_bin_limit) {
        vector<vector<double>> overlap_vec(num_overlaps, vector<double>(2));
        for (int i = 0; i < num_overlaps; i++) {//初始化
            overlap_vec[i][0] = err_sum[i];
            overlap_vec[i][1] = overlap_lines[i];
        }
        sort(overlap_vec.begin(), overlap_vec.end());

        for (auto iter_overlaps = overlap_vec.begin(); num_bins < min_bin_limit; num_bins++, iter_overlaps++) {
            *iter_lines = (*iter_overlaps)[1];
        }
        num_bins = min_bin_limit;
//        auto posi = min_element(err_sum.begin(), err_sum.end());
//        lines[0] = (overlap_lines[posi - err_sum.begin()]);
//        //cout << "Final line is: " << lines[0] << endl;
    }
    lines.resize(num_bins - 1);
    return lines;
}

vector<vector<double> > prep_index(vector<double> &c1, vector<double> &c2, double median_1, double median_2) {
    // sort two clusters
    sort(c1.begin(), c1.end());
    sort(c2.begin(), c2.end());
    auto c1_l = lower_bound(c1.begin(), c1.end(), median_1);
    auto c1_r = upper_bound(c1.begin(), c1.end(), median_2);
    vector<double> c1_new(c1_l, c1_r);// part of clsuter 1 in between median_1 and median_2

    auto c2_l = lower_bound(c2.begin(), c2.end(), median_1);
    auto c2_r = upper_bound(c2.begin(), c2.end(), median_2);
    vector<double> c2_new(c2_l, c2_r);// part of clsuter 2 in between median_1 and median_2

    if (c1_new.size() + c2_new.size() == 0) {
        vector<vector<double> > mean_index = vector<vector<double> >(3, vector<double>(0, 0));
        return mean_index;
    }

    // merge them together
    auto iter_c1 = c1_new.begin();
    auto iter_c2 = c2_new.begin();
    vector<vector<double> > data_index = vector<vector<double> >(3, vector<double>(c1_new.size() + c2_new.size(), 0));

    auto real_data = data_index[0].begin();
    auto c1_index_iter = data_index[1].begin();
    auto c2_index_iter = data_index[2].begin();

    int prev_1 = distance(c1_l, c1.end()) + 1;
    int prev_2 = distance(c2.begin(), c2_l);
    int final_index_c1 = distance(c1_r, c1.end());
    int final_index_c2 = distance(c2.begin(), c2_r - 1);

    for (; real_data != data_index[0].end(); ++real_data, ++c1_index_iter, ++c2_index_iter) {
        if (iter_c1 == c1_new.end()) {//c1 is empty, put c2
            *real_data = *iter_c2;
            iter_c2++;
            prev_2++;
            prev_1 = final_index_c1;

        } else if (iter_c2 == c2_new.end()) {//c2 is empty, put c1
            *real_data = *iter_c1;
            prev_1--;
            iter_c1++;
            prev_2 = final_index_c2;

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

        *c1_index_iter = prev_1;
        *c2_index_iter = prev_2;
    }

    vector<vector<double> > mean_index = vector<vector<double> >(3, vector<double>(data_index[0].size() - 1, 0));

    for (int i = 0; i < int(data_index[0].size()) - 1; ++i) {
        double mean = (data_index[0][i] + data_index[0][i + 1]) / 2.0;
        int c1_index = data_index[1][i + 1];
        int c2_index = data_index[2][i];
        mean_index[0][i] = mean;
        mean_index[1][i] = c1_index;
        mean_index[2][i] = c2_index;
    }
    return mean_index;
}

double binary_search_index(const vector<vector<double> > &c_index, const int left, const int right, const int size_c1,
                           const int size_c2, bool &overlap, vector<double> &err_sum) {

    if (right < left) {
        throw "VALUE ERROR!\nRight bound less than left bound!\n";
    }

    int line_index = ceil(float(left + (right - left) / 2.0));
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
