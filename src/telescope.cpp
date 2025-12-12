//
// Created by Jiandong Wang on 2023/4/29.
//

#include "telescope.h"

using namespace std;

static const vector<double> *spx_c;

static bool compi_c(size_t i, size_t j) {
   return (*spx_c)[i] < (*spx_c)[j];
}

vector<size_t> tele(const vector<vector<double>> &data, const size_t n_lines) {
    size_t n_dim = data.size();
    size_t n_points = data[0].size();

    vector<vector<double>> lines(n_dim, vector<double>(n_lines, 0));
    vector<vector<size_t>> tags_num_by_dim(n_dim, vector<size_t>(n_points, 0));
    vector<vector<size_t>> tags_num(n_points, vector<size_t>(n_dim, 0));

    split(data, n_lines, lines, tags_num_by_dim, tags_num);

    HashTableForTele htt(tags_num);


    //htt.print();

    auto labels = htt.gen_label();


    //vector<size_t> labels = gen_label(hct, n_points);
    return labels;
    //todo: traversal the whole table by key and assign the label correspondingly
}


void split(const vector<vector<double>> &data, const size_t n_lines, vector<vector<double>> &lines,
           vector<vector<size_t>> &tags_num_by_dim, vector<vector<size_t>> &tags_num) {
    size_t n_dim = data.size();
    size_t n_points = data[0].size();
    double scope = n_points / (double) (n_lines + 1);

    //vector<string> tags(n_points, " ");
    /**
     * @note modified by J.Wang 2023.05.07 \n
     *       change the user input from scope (number of points in each zone) to number of lines.
     *       cutting by #lines is more flexible and #points in each zone will be more balanced (the difference no larger than 1)
     */
    /*for (size_t d = 0; d < n_dim; ++d) {
        // loop on dimension
        size_t index = 0;
        pre_point = 0;
//        auto iter_lines = lines[d].begin();
        size_t label = 0;
        auto iter_tag = tags_num[d].begin();
        for (auto x: data[d]) {
            ++index;
            if (index == pre) {
                pre_point = x;
            } else if (index == scope) {
//                *iter_lines = (pre_point + x) / 2;
//                ++iter_lines;

                lines[d][label] = (pre_point + x) / 2;
                ++label;
                index = 0;
            }
            *iter_tag = label;
            ++iter_tag;
        }
    }*/
    /**
     * @note modified by J.Wang 2023.05.07 \n
     *       only calculate points on the left and right side of each line instead of looping the whole dataset to save the running time
     *  todo: test the double division, check +1 or -1, and the range
     *  todo: maybe the line should start from 1 instead of 0;
     */

    auto order = sort_data(data);
    auto lines_iter = lines.begin();
    for (size_t d = 0; d < n_dim; ++d) {
        //for each dimension
/*
 *      commented, goes wrong when the last scope is not full(e.g. 11 point and scope is 3, the last one only has 2 points)for (int i = 1; i <= n_lines; ++i) { // calculate line values
 *
            int index_r = scope * i; // floor cast, e.g. 2.7 -> 2
            double line_val = (data[d][index_r - 1] + data[d][index_r]) / 2;
            (*lines_iter)[i] = line_val;
            // assign tag to tag_num
            vector<size_t> tags(scope, i);
            tags_num[i] = tags;
        }*/
        // calculation for line position
        int index = 0;
        int right_side = ceil((index + 1) * scope);
        int left_side = right_side - 1;

        for (size_t i = 0; i < n_points; ++i) { // calculate line values
            if (i == (size_t) right_side) {
                double line_val = (data[d][order[d][left_side]] + data[d][order[d][right_side]]) / 2;
                (*lines_iter)[index] = line_val;
                ++index;
                right_side = ceil((index + 1) * scope);
                left_side = right_side - 1;
            }

            //assign tag to each point
            tags_num[order[d][i]][d] = index;
            tags_num_by_dim[d][order[d][i]] = index;
        }
        ++lines_iter;
    }
    //return tags_num;
    //todo: Process the tag_num, transform it to tag, string.
    // Maybe do this in another function
    // Maybe can using vector<int> as key directly, instead of transforming into string
}

vector<vector<size_t>> split(const vector<vector<double>> &data, const size_t n_lines) {
    /**
     * @param data the dataset that will be processed, organized as \<dims \<points\>\>
     * @param scope an integer indicate how many points are allowed in each zone
     */
    size_t n_dim = data.size();
    // size_t n_points = data[0].size();
    //size_t n_lines = n_points / scope;
    //size_t pre = scope - 1;
    // double scope = n_points / n_lines;

    //double pre_point;
    vector<vector<double>> lines(n_dim, vector<double>(n_lines, 0));
    vector<vector<size_t>> tags_num_by_dim(n_dim, vector<size_t>(n_lines, 0));
    vector<vector<size_t>> tags_num(n_lines, vector<size_t>(n_dim, 0));

    split(data, n_lines, lines, tags_num_by_dim, tags_num);
    return tags_num;
    //todo: Process the tag_num, transform it to tag, string.
    // Maybe do this in another function
    // Maybe can using vector<int> as key directly, instead of transforming into string
}

/*vector<size_t> gen_label(const HCountTable hct, size_t n) {
    vector<size_t> labels(n, 0);
    int i = 0;
    for (auto & p : hct.m_table) {
        for (auto x = 0; x < p.second.posi; ++x) {
            labels[p.second.by_index[x]] = i;
        }
        ++i;
    }
}*/

vector<vector<double>> transform_by_dim(const vector<vector<double>> &data) {
    /**
     * @brief transform the data from \<points \<dims\>\> to \<dims \<points\>\>
     *
     * @param data the dataset that organized as \<points \<dims\>\>
     *
     * @return the transformed dataset that organized as \<dims \<points\>\>
     */
    size_t n_points = data.size();
    size_t n_dim = data[0].size();
    vector<vector<double>> data_by_dim(n_dim, vector<double>(n_points, 0));
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t d = 0; d < n_dim; ++d) {
            data_by_dim[d][i] = data[i][d];
        }
    }
    return data_by_dim;
}

vector<vector<size_t>> sort_data(const vector<vector<double>> &data) {
    size_t n_dims = data.size();
    size_t n_points = data[0].size();

    auto order = vector<vector<size_t>>(n_dims, vector<size_t>(n_points, 0));
    vector<size_t> order_sample(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        order_sample[i] = i;
    }
    for (size_t dim = 0; dim < n_dims; ++dim) {
        spx_c = &data[dim];
        auto data_order = order_sample;
        sort(data_order.begin(), data_order.end(), compi_c);
        order[dim] = data_order;
    }
    return order;
}
//
//void gen_tag(const vector<string> &tag, const vector<vector<size_t>> &tag_num) {
//    /**
//     * @param tag a vector of string with the size of dataset, each element is a label(tag) for the point, will be used as the key in hashtable
//     * @param tag_num a numerical matrix organized as \<dims \<points\>\>, indicated which zone the point belongs to
//     */
//    for (auto tag_d: tag_num) {
//        auto iter_tag_n = tag_d.begin();
//        auto iter_tag = tag.begin();
//        for (auto tag_d_p: tag_d) {
//            *iter_tag = strcpy(*iter_tag, to_string(tag_d_p));
//        }
//    }
//}
