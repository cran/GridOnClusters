//
// Created by Jiandong Wang on 2022/1/16.
//

#include <R.h>
#include "Cutting_Cluster_dp.h"

/*vector<int> Cutting_Cluster_1D(vector<int> &X, unsigned int K, unsigned int Q_max, unsigned int Q_min) {
    /// min_bin_limit is the Q_min
//    int n = X.size();

    vector<int> posi = Cutting_Cluster(X, K, Q_max, Q_min);

   //calculate the index
    vector<vector<int>> index(K, vector<int>(n, 0));
    index = prep_index(X, K);

    vector<double> lines(posi.size() - 1);

    for (int i = 0; i < posi.size() - 1; ++i) {
        int posi_a = posi[i];
        int posi_b = posi_a + 1;

        int label_a = X[posi_a];
        int label_b = X[posi_b];

        int index_a = index[label_a - 1][posi_a];
        int index_b = index[label_b - 1][posi_b];

        double val_a = cluster.get_points(label_a, dim_input, index_a);
        double val_b = cluster.get_points(label_b, dim_input, index_b);

        double mid = (val_a + val_b) / 2.0;
        lines[i] = mid;
    }

    return posi;
}*/


vector<int> Cutting_Cluster(vector<int> &X, unsigned int K, unsigned int Q_max, unsigned int Q_min, double cut_off) {
    /**
     * @param X is a vector of cluster label\n
     * @param k is the number of clusters\n
     * @param Q_min indicate the lower bound of Q\n
     * @param Q_max indicate the upper bound of Q, should be greater or equals to Q_min.
     */


    vector<int> init(K, 0);
    size_t n = X.size();
    vector<vector<double>> S(n, vector<double>(Q_max, 0)); // score matrix
    vector<vector<int>> P(n, vector<int>(Q_max, 0)); // path matrix

    /// initialize the first base case
    for (size_t i = 0; i < n; ++i) {
        init[X[i]] += 1;
        S[i][0] = *max_element(init.begin(), init.end());
    }

    ///calculate the index
    vector<vector<int>> index(K, vector<int>(n, 0));
    index = prep_index(X, (int) K);


    /// calculate score
    for (size_t q = 1; q < Q_max; ++q) { // loop for col
        for (size_t i = q; i < n; ++i) { // loop for row
            double max_score = 0;
            int pre_node = 0;
            double temp_score;

            for (size_t j = q - 1; j < i; ++j) { // loop for max score select
                unsigned int max_score_region = 0;
                for (size_t x = 0; x < K; ++x) {
                    unsigned int temp_s = index[x][i] - index[x][j];
                    if (max_score_region < temp_s) {
                        max_score_region = temp_s;
                    }
                }
                temp_score = S[j][q - 1] + max_score_region;

                if (max_score < temp_score) {
                    max_score = temp_score;
                    pre_node = (int) j;// pre is the end of last cluster
                }
            }
            S[i][q] = max_score;
            P[i][q] = pre_node;
        }
    }

    vector<int> best_path;
    vector<vector<int>> best_table;
    //mydouble best_estimate = 0;
    //mydouble best_Fstat;
    //double best_pval;
    double max_upsilon = 0;

    vector<vector<int>> all_path(Q_max - Q_min + 1, vector<int>(Q_max, -1));
    vector<double> all_upsilon(Q_max - Q_min + 1, -1);

    for (unsigned int q = Q_min; q <= Q_max; ++q) {
        vector<int> path;
        vector<vector<int>> table;
        double upsilon;

        path = trace_back(q - 1, P);
        table = prep_table(X, path, Q_max, K);

        upsilon = upsilon_stat(table);
        if (upsilon > max_upsilon) {
            best_path = path;
            max_upsilon = upsilon;
            best_table = table;
        }
        all_path[q - Q_min] = path;
        all_upsilon[q - Q_min] = upsilon;
    }
    int df = (best_table[0].size() - 1) * (best_table.size() - 1);
    //double median_pdf = df * (1 - pow((2.0 / (9.0 * df)), 3));
    auto pval = (double) ChisqPvalue(max_upsilon, df);

    /// if not satisfy the customer require
    if (pval > cut_off) {
        best_path = all_path[0];
        //max_upsilon = all_upsilon[0];
    }

/*    if (median_pdf > max_upsilon) {
        best_path = all_path[0];
        max_upsilon = all_upsilon[0];
    }*/

    return best_path;
}

vector<int> trace_back(unsigned int start_k, const vector<vector<int>> &P) {
    vector<int> path(start_k + 1, 0);
    path[start_k] = (int) P.size() - 1;
    if (start_k == 0) {
        return path;
    }
    path[start_k - 1] = P[P.size() - 1][start_k];
    for (size_t i = start_k - 1; i > 0; --i) {
        path[i - 1] = P[path[i]][i];
    }
    return path;
}

vector<vector<int>> prep_table(const vector<int> &X, const vector<int> &stop_posi, unsigned int Q, unsigned int K) {
    // funchi table, row is cluster label, col is q
    int max_K = max((int)Q, (int)K);
    vector<vector<int>> table(K, vector<int>(max_K, 0));
    auto posi_iter = stop_posi.begin();
    int q = 0;
    for (size_t i = 0; i < X.size(); ++i) {
        ++table[X[i]][q];
        if (i == (size_t)*posi_iter) {
            ++posi_iter;
            ++q;
        }
    }
    return table;
}


void print_table(const vector<vector<int>> &table) {
    size_t nrow = table.size();
    size_t ncol = table[0].size();

    Rprintf("\t");     // cout << "\t";
    
    for (size_t j = 0; j < ncol; ++j) {
      Rprintf("c%zu\t", j);
      // cout << "c" << j << "\t";
    }

    for (size_t i = 0; i < nrow; ++i) {
        // cout << endl << i << "\t";
        Rprintf("\n%zu\t", i);
      
        for (size_t j = 0; j < ncol; ++j) {
            // cout << table[i][j] << "\t";
            Rprintf("%d\t", table[i][j]);
        }
    }
    Rprintf("\n"); //cout << endl;
}

vector<vector<int>> prep_index(vector<int> &X, int K) {
    size_t n = X.size();
    vector<vector<int>> index(K, vector<int>(n, 0));

    for (size_t i = 0; i < n; ++i) {
        //calculate the index
        for (int j = 0; j < K; ++j) {
            if (i > 0) {
                index[j][i] = index[j][i - 1];
            } else {
                index[j][i] = 0;
            }
        }
        index[X[i]][i] += 1;
    }
    return index;
}
