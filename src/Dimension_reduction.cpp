//
// Created by Jiandong Wang on 2022/11/3.
//

#include "Dimension_reduction.h"
#include <iomanip>

vector<bool> dim_reduct(
    Cluster &clusters, 
    vector<vector<double>> & lines_removed, 
    vector<vector<size_t>> &best_table,
    const double global_BIC, 
    double &best_BIC, size_t n_thread, size_t step) {
    /**
     * The main function to do dimension reduction for the raw data.
     * @param lines_removed separate lines for each dimension after the dimension reduction.
     * @param best_table the best contig-table for data after dimension reduction.
     * @param global_BIC the best BIC for original data, the smaller the better.
     * @return a vector of bool indicate which dimension has been removed
     */

    // size_t N = clusters.get_sample_size();
    // size_t K = clusters.get_K();
    size_t n_dim = clusters.get_dims();

    double global_best_BIC = global_BIC;
    vector<bool> dimensions_removed(n_dim, false);
    vector<size_t> num_lines(n_dim, 0);
    //vector<vector<size_t>> best_table;

    for (size_t i = 0u; i < n_dim; ++i) {
        num_lines[i] = lines_removed[i].size();
    }
    //auto lines_removed = lines_all_dim;

    for (auto i = 0u; i < n_dim; ++i) {
        //out loop for dimension

        vector<double> BIC_collection(n_dim, numeric_limits<double>::lowest());
        vector<vector<size_t>> local_table;
        // table with the worst(largest) BIC score
        double local_best_BIC = numeric_limits<double>::max();
        size_t best_index = 0;

        for (size_t j = 0; j < n_dim; ++j) {
            // inner loo[ for dimension
            // calculate the BIC for all cases of j-th dim removed
            if (dimensions_removed[j]) {
                //this dimension has already been removed
                continue;
            }
            vector<vector<size_t>> conti_table;
            double bic_temp = dim_BIC(clusters, lines_removed, num_lines, j, conti_table, n_thread, step);
            BIC_collection[j] = bic_temp;
            if (local_best_BIC > bic_temp) {
                //update the local worst(smallest) BIC
                local_best_BIC = bic_temp;
                //BIC_collection[j] = bic_temp;
                local_table = conti_table;
                best_index = j;
            }
        }

        //size_t removed_index = max_element(BIC_collection.begin(), BIC_collection.end()) - BIC_collection.begin();
        if (global_best_BIC > BIC_collection[best_index]) {
            // if the BIC of the data that without this dimension is better than the global one, which means the data
            // would be better by removing this dimension
            global_best_BIC = local_best_BIC;
            // remove that dimension
            dimensions_removed[best_index] = true;
            num_lines[best_index] = 0;
            lines_removed[best_index].clear();
            best_table = local_table;
            best_BIC = BIC_collection[best_index];
        }
        //cout << endl;
    }
    return dimensions_removed;
}

double
dim_BIC(Cluster &clusters, vector<vector<double>> lines_dims, vector<size_t> num_lines, size_t dim_removed,
        vector<vector<size_t>> &conti_table, size_t n_thread, size_t step) {
    /**
     * It will take the raw data and the lines to calculate how many points has fallen in between two lines for all
     * dimensions. Generate table accordingly, and calculate the BIC for that table.
     */

    size_t N = clusters.get_sample_size();

    //remove all lines on the dimension that need be removed
    lines_dims[dim_removed].clear();
    num_lines[dim_removed] = 0;

    //generate discretization table
    vector<vector<int>> discr_data = discretize_data(clusters, lines_dims, n_thread, step, false);

    //hash discr_table into hash table
    vector<size_t> labels_t(N, 0);
    auto iter_lab = labels_t.begin();
    for (int label: clusters.get_labels()) {
        *iter_lab = abs(label);
        ++iter_lab;
    }

    MyHash hct(discr_data, labels_t);

    //generate conti-table
    //Conti-table: row is the cluster(predict), col is the class(truth),

    conti_table = hct.gen_table();

    //run BIC-for-cluster for dim-reduction

/*
    // convert to transform table
    vector<vector<size_t>> conti_table_t(conti_table[0].size(), vector<size_t>(conti_table.size(), 0));

    for (int i = 0; i < conti_table.size(); ++i) {
        for (int j = 0; j < conti_table[0].size(); ++j) {
            conti_table_t[j][i] = conti_table[i][j];
        }
    }
*/
    double BIC_table = BIC_for_table(conti_table, N, num_lines);
    return BIC_table;
}

double BIC_for_table(vector<vector<size_t>> table, size_t N, vector<size_t> nlines) {
    /// @param table is a two dimension q x k vector with cluster as row and class as col.
    size_t Q = table.size();
    size_t K = table[0].size();
    size_t n_cubes = table.size();
    // size_t totle_cubes = 1;
    // size_t emp_cubes = 0;

    // size_t total_for_col = 0;
    double LL = 0;

/*    for (size_t k = 0; k < K; ++k) {
        /// loop for classes
        for (size_t q = 0; q < Q; ++q) {
            /// loop for clusters
            if (table[q][k] == 0) {
                continue;
            }
            total_for_col += table[q][k];
            LL += table[q][k] * log(table[q][k]);
        }
    }*/
    for (size_t q = 0; q < Q; ++q) {
        /// loop for clusters
        size_t sum_cube = 0;
        for (size_t k = 0; k < K; ++k) {
            /// loop for classes
            if (table[q][k] == 0) {
                continue;
            }
            // total_for_col += table[q][k];
            LL += table[q][k] * log(table[q][k]);
            sum_cube += table[q][k];
        }
        LL = LL - sum_cube * log(sum_cube);
    }

    size_t sum_num_lines = 0;
    for (size_t nline: nlines) {
        sum_num_lines += nline;
        // totle_cubes *= (nline + 1);
    }
    // emp_cubes = totle_cubes - n_cubes;

    //double log_likelihood = LL - LL * log(LL);
    double log_likelihood = LL;
    //double BIC = -2 * log_likelihood + (Q - 1 + (K - 1) * Q) * log((double) N);
    //double BIC = -2 * log_likelihood + n_cubes * log((double) N);
    double BIC = -2 * log_likelihood + (sum_num_lines + n_cubes * (K - 1)) * log((double) N);

/*    //print out
    for (size_t nline: nlines) {
        cout << nline << " ";
    }
    cout << setw(14) << log_likelihood << setw(13) << sum_num_lines << setw(13) << totle_cubes << setw(13) << n_cubes
         << setw(13) << BIC << ";"
         << endl;*/

    return BIC;
}
