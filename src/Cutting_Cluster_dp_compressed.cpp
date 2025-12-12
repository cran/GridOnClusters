/**
 * @file Cutting_Cluster_dp_compressed.cpp
 * @author Jiandong Wang
 * @date 3/08/22
 * @copyright NMSU Song lab
 */

#include "Cutting_Cluster_dp_compressed.h"

vector<int> Cutting_Cluster_dp_compressed(vector<int>& labels, vector<double>& data, 
                                          size_t K, size_t Q_max, size_t Q_min, 
                                          double & best_bic, const string& method, 
                                          const bool entropy) {
    /**
     *  @param labels the labels with corresponding order of the data in that dimension\n
     *  @param data the data points in the increasing order of that dimension, store the data by cluster\n
     *  @param K the number of the clusters\n
     *  @param best_BIC the best BIC score for this dimension
     *  @param Q_max and \b Q_min indicate the upper and lower bound of the \b Q, where \b Q is the number of bins in this dimension\n
     *  @param method the discretization method selected\n
     *  @param entropy whether using entropy or likelihood\n
     */
    size_t N = labels.size();
    vector<int> labels_compressed(N), labels_weights(N), weight_sum(N);
    vector<vector<int>> index(K, vector<int>(N, 0));

    // Compress input labels into blocks for faster processing
    compress_labels(labels, labels_compressed, labels_weights, weight_sum, index, N, K);

    /// dynamic programming initialize
    int n_blocks = static_cast<int>(labels_compressed.size());
    Q_max = min(static_cast<int>(Q_max), n_blocks); // Clamp Q_max to the number of available blocks

    // Initialize dynamic programming matrices
    vector<vector<double>> S(n_blocks, vector<double>(Q_max, 0.0)); // Score matrix
    vector<vector<int>> P(n_blocks, vector<int>(Q_max, 0)); // Path matrix
    vector<vector<double>> S_e(n_blocks, vector<double>(Q_max, 0.0)); // Entropy score matrix

    // Score computation based on selected method
    if (method == "DP exact likelihood") {
        cal_score_binomial(data, labels_compressed, labels_weights, S, S_e, P, index, Q_max, K, n_blocks, entropy);
    } else if (method == "DP approx likelihood"|| method == "DP approx likelihood two way") {
        cal_score_binomial_recur_main(data, labels_compressed, labels_weights, S, P, index, Q_max, K, n_blocks,
                                      entropy);
    } else if (method == "DP compressed majority") {
        cal_score_absolute(labels_compressed, labels_weights, S, P, index, Q_max, K, n_blocks);
    }

    //print_out_table(S);
    vector<int> best_path;
    vector<vector<int>> best_table;
    best_bic = numeric_limits<double>::infinity();

    vector<vector<int>> all_path_real(Q_max - Q_min + 1, vector<int>(Q_max, -1));
    vector<double> all_BIC(Q_max - Q_min + 1, numeric_limits<double>::infinity());

    for (unsigned int q = Q_min; q <= Q_max; ++q) {
        vector<int> comp_path = trace_back_dp(q - 1, P);
        vector<int> real_path(q);
        vector<vector<int>> table;
        double bic;

        //path matrix store the index of the end of clusters (index starts with 0)
        comp_path = trace_back_dp(q - 1, P);

        // Convert compressed path to real data positions
        vector<vector<int>> index_over_channel(K, vector<int>(comp_path.size(), 0));
        for (size_t i = 0; i < comp_path.size(); ++i) {
            int path_index = comp_path[i];
            real_path[i] = weight_sum[path_index] - 1;
            for (size_t j = 0; j < K; ++j) {
                index_over_channel[j][i] = index[j][path_index];
            }
        }

        table = prep_table_dp(labels_compressed, labels_weights, comp_path, q, K);

        if (method == "DP exact likelihood" || method == "DP approx likelihood"|| method == "DP approx likelihood two way") {
            bic = entropy ? evaluation_entropy(N, q, K, S, S_e) : evaluation_likelihood(N, q, K, S);
        } else if (method == "DP compressed majority") {
            bic = evaluation_counts(data, labels, index_over_channel, real_path, table);
        } else {
            bic = 0.0;
        }


        // Select the best path based on BIC
        if (bic < best_bic) {
            best_path = real_path;
            best_bic = bic;
            best_table = table;
        }
        all_BIC[q - Q_min] = bic;
        all_path_real[q - Q_min] = real_path;
    }

    //print_out_table(S);
    //print_out_table(P);
    //
    // cout << "best BIC is :  " << best_bic << endl;
    return best_path;
}

/*
void compress_labels(const vector<int>& labels, vector<int>& labels_compressed, vector<int>& labels_weights,
                     vector<int>& weight_sum, vector<vector<int>>& index, unsigned int N, unsigned int K)
{
    /// compress the label
    int cur_label = labels[0];
    int cur_weight = 1;
    int cur_index = 0;

    for (size_t i = 1; i < N; ++i)
    {
        if (cur_label == labels[i])
        {
            cur_weight++;
        }
        else
        {
            // label changed, no longer in the same cluster, update the variable
            labels_compressed[cur_index] = cur_label;
            labels_weights[cur_index] = cur_weight;
            weight_sum[cur_index] = (int)i;
            update_index(index, cur_index, cur_weight, cur_label, K);
            cur_label = labels[i];
            cur_weight = 1;
            cur_index++;
        }
    }
    {
        labels_compressed[cur_index] = cur_label;
        labels_weights[cur_index] = cur_weight;
        weight_sum[cur_index] = weight_sum[cur_index - 1] + cur_weight;
        update_index(index, cur_index, cur_weight, cur_label, K);
        cur_index++;
    }
    labels_weights.resize(cur_index);
    labels_compressed.resize(cur_index);
    weight_sum.resize(cur_index);
    for (auto iter = index.begin(); iter != index.end(); ++iter)
    {
        (*iter).resize(cur_index);
    }
}
*/

void compress_labels(const vector<int>& labels, vector<int>& labels_compressed,
                     vector<int>& labels_weights, vector<int>& weight_sum,
                     vector<vector<int>>& index, const unsigned int N, const unsigned int K) {
    ///
    /// @param N the number of the points\n
    /// @param K the number of the clusters\n
    ///
    if (labels.empty()) return; // Handle empty input case

    labels_compressed.resize(N); // Preallocate max possible size
    labels_weights.resize(N);
    weight_sum.resize(N);

    int cur_label = labels[0];
    int cur_weight = 1;
    int cur_index = 0;

    for (size_t i = 1; i < N; ++i) {
        if (labels[i] == cur_label) {
            cur_weight++;
        } else {
            // Store compressed label information
            labels_compressed[cur_index] = cur_label;
            labels_weights[cur_index] = cur_weight;
            weight_sum[cur_index] = i; // Current index as weight sum

            update_index(index, cur_index, cur_weight, cur_label, K);

            // Reset for new label
            cur_label = labels[i];
            cur_weight = 1;
            cur_index++;
        }
    }

    // Store the last group
    labels_compressed[cur_index] = cur_label;
    labels_weights[cur_index] = cur_weight;
    weight_sum[cur_index] = (cur_index > 0) ? weight_sum[cur_index - 1] + cur_weight : cur_weight;
    update_index(index, cur_index, cur_weight, cur_label, K);

    // Trim vectors to actual size
    labels_compressed.resize(cur_index + 1);
    labels_weights.resize(cur_index + 1);
    weight_sum.resize(cur_index + 1);

    // Resize index matrix
    for (auto& row : index) {
        row.resize(cur_index + 1);
    }
}

void update_index(vector<vector<int>>& index, int cur_index, int cur_weight, int cur_label, unsigned int K) {
    // Copy previous counts if not the first block
    if (cur_index > 0) {
        for (size_t i = 0; i < K; ++i) {
            index[i][cur_index] = index[i][cur_index - 1];
        }
    } else {
        // First block: initialize all counts to 0
        for (size_t i = 0; i < K; ++i) {
            index[i][cur_index] = 0;
        }
    }

    // Increment count for the current label
    index[cur_label][cur_index] += cur_weight;
}

vector<int> trace_back_dp(unsigned int start_k, const vector<vector<int>>& P) {
    // Initialize the path vector with size start_k + 1
    vector<int> path(start_k + 1);

    // Start from the last bin, which ends at the final row of P
    path[start_k] = static_cast<int>(P.size()) - 1;

    // Trace back through the path matrix
    for (int i = static_cast<int>(start_k); i > 0; --i) {
        path[i - 1] = P[path[i]][i];
    }
    return path;
}

vector<vector<int>> prep_table_dp(const vector<int>& label, const vector<int>& weight, const vector<int>& stop_posi,
                                  unsigned int Q,
                                  unsigned int K) {
    // funchi table, row is cluster label, col is q
    //Padding zero to have a fixed size table (degrees of freedom), for upsilon test, etc. (input using Q_max)
    //int max_K = max(Q, K);

    // Create a K x Q zero-initialized matrix
    vector<vector<int>> table(K, vector<int>(Q, 0));

    // Iterator over stopping positions
    auto stop_it = stop_posi.begin();
    unsigned int current_bin = 0;

    for (size_t i = 0; i < label.size(); ++i) {
        // Accumulate weight for each cluster in the current bin
        table[label[i]][current_bin] += weight[i];

        // Move to next bin if current index matches a stopping position
        if (stop_it != stop_posi.end() && i == static_cast<size_t>(*stop_it)) {
            ++stop_it;
            ++current_bin;
            if (current_bin >= Q) break; // Avoid overflow in case of malformed input
        }
    }
    return table;
}

void cal_score_absolute(const vector<int>& labels_compressed, const vector<int>& labels_weights,
                        vector<vector<double>>& S,
                        vector<vector<int>>& P, const vector<vector<int>>& index, unsigned int Q_max, unsigned int K,
                        int n_blocks) {
    /// initialize the first base case\n
    ///\b m is the number of the compressed regions
    vector<int> init(K, 0);
    // maybe can initialize index here to avoid unnecessary traversal (resize)
    for (int i = 0; i < n_blocks; ++i) {
        init[labels_compressed[i]] += labels_weights[i];
        int max_init = *max_element(init.begin(), init.end());
        S[i][0] = max_init;
    }

    for (size_t q = 1; q < Q_max; ++q) {
        // loop for col
        for (size_t i = q; i < (size_t)n_blocks; ++i) {
            // loop for row
            double max_score = numeric_limits<double>::min();
            int pre_node = 0;
            double temp_score;

            for (size_t j = q - 1; j < i; ++j) {
                // loop for max score select
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
                    pre_node = j; // pre is the end of last cluster
                }
            }
            S[i][q] = max_score;
            P[i][q] = pre_node;
        }
    }
}

void cal_score_binomial(const vector<double>& data, const vector<int>& labels_compressed,
                        const vector<int>& labels_weights, vector<vector<double>>& S, vector<vector<double>>& S_e,
                        vector<vector<int>>& P, const vector<vector<int>>& index, unsigned int Q_max, unsigned int K,
                        int n_blocks, bool entropy) {
    // Initialize the base case for the first bin
    vector<int> init_weights(K, 0);
    vector<int> total_weights(n_blocks, 0);

    // maybe can initialize index here to avoid unnecessary traversal (resize)
    for (int i = 0; i < n_blocks; ++i) {
        init_weights[labels_compressed[i]] += labels_weights[i];

        // Compute cumulative weights
        total_weights[i] = (i == 0) ? labels_weights[i] : total_weights[i - 1] + labels_weights[i];

        double score = 0;
        // Compute log-likelihood score for each cluster
        for (size_t x = 0; x < K; ++x) {
            if (init_weights[x] == 0) continue;
            score += init_weights[x] * log((double)init_weights[x] / total_weights[i]);
        }

        S[i][0] = score;
        P[i][0] = i;
        if (entropy) {
            S_e[i][0] = score / total_weights[i];
        }
    }

    // Dynamic programming to fill out remaining bins
    for (size_t q = 1; q < Q_max; ++q) {
        // loop for col, q is the bin to be filled
        for (size_t i = q; i < (size_t)n_blocks; ++i) {
            // loop for row, block i to be filled
            double max_score = -1 * numeric_limits<double>::infinity(); // DBL_MIN;
            double max_score_etropy = -1 * numeric_limits<double>::infinity(); // DBL_MIN;
            int pre_node = 0;

            for (size_t j = q - 1; j <= i; ++j) {
                // loop for max score select, j is the start of that region
                if (j > 0 && j != q - 1 && data[j] == data[j - 1]) continue;
                // Skip redundant splits
                // In high-dimensional data, it’s possible for data[j] and data[j-1] to have identical values yet
                // belong to different clusters. Placing a cut here would be misleading and ambiguous.

                /*
                unsigned int total_points = 0;
                double score_region = 0;

                // Calculate region likelihood score
                for (size_t x = 0; x < K; ++x) {
                    // loop for all cluster in that region
                    unsigned int count = index[x][i] - index[x][j];
                    if (count == 0) { continue; }
                    total_points += count;
                    score_region += count * log(count);
                }
                score_region = score_region - total_points * log(total_points);
                */


                double cll = loglikelihood(j + 1, i, index, entropy);
                double temp_score = S[j][q - 1] + cll;

                if (max_score < temp_score) {
                    max_score = temp_score;
                    pre_node = j; // pre is the end of last cluster
                }

                // Update entropy score if enabled
                if (entropy) {
                    //score_region = score_region / total_points;
                    temp_score = S_e[j][q - 1] + cll;
                    if (max_score_etropy < temp_score) {
                        max_score_etropy = temp_score;
                        pre_node = j; // pre is the end of last cluster
                    }
                }
            }

            // Store best scores and backtrack pointers
            S[i][q] = max_score;
            P[i][q] = pre_node; // path store the beginning of the cluster     8/2/2022
            if (entropy) {
                S_e[i][q] = max_score_etropy;
            }
        }
    }
}

void cal_score_binomial_recur_main(const vector<double>& data,
                                   const vector<int>& labels_compressed,
                                   const vector<int>& labels_weights,
                                   vector<vector<double>>& S,
                                   vector<vector<int>>& P,
                                   const vector<vector<int>>& index,
                                   unsigned int Q_max,
                                   unsigned int K,
                                   int n_blocks,
                                   bool entropy) {
    // Base case q = 0
    vector<int> init_weights(K, 0);
    vector<int> total_weights(n_blocks, 0);

    for (int i = 0; i < n_blocks; ++i) {
        init_weights[labels_compressed[i]] += labels_weights[i];
        total_weights[i] = (i == 0) ? labels_weights[i] : total_weights[i - 1] + labels_weights[i];

        /*double score = 0;
        for (size_t x = 0; x < K; ++x) {
            if (init_weights[x] > 0)
                score += init_weights[x] * log(static_cast<double>(init_weights[x]) / total_weights[i]);
        }
        S[i][0] = score;
        */
        double score = 0.0L;
        for (size_t x = 0; x < K; ++x) {
            if (init_weights[x] > 0)
                score += init_weights[x] * log(double(init_weights[x]) / total_weights[i]);
            // cout << std::fixed << std::setprecision(20) << "init_weights:\t" << init_weights[x] << "\t total_weights:\t"
            //     << total_weights[i] << "\t frac:\t" << double(init_weights[x]) / total_weights[i] << "\t log(frac):\t" <<
            //         log(double(init_weights[x]) / total_weights[i]) << "\t score:\t" << score << endl;
        }
        S[i][0] = static_cast<double>(score);

        P[i][0] = i;
    }

    // DP for q > 0
    for (size_t q = 1; q < Q_max; ++q) {
        // Only start from i = q (need at least q+1 points for q+1 bins)
        if (n_blocks > (int) q)
            cal_score_binomial_recur(q, static_cast<int>(q), 0, n_blocks - 1, n_blocks - 1, S, P, index, data, entropy);
    }
}

void cal_score_binomial_recur(int q,
                              int il, int jl,
                              int iu, int ju,
                              vector<vector<double>>& S,
                              vector<vector<int>>& P,
                              const vector<vector<int>>& index,
                              const vector<double>& data,
                              bool entropy) {
    if (il > iu) return;

    int i = (il + iu) / 2;
    if (i < q) {
        // Not enough points to fill q+1 bins
        S[i][q] = 0.0; // or -inf if you prefer invalid
        P[i][q] = -1;
        return;
    }

    if (jl < q - 1) {
        jl = q - 1;
    }

    double best_score = -1 * numeric_limits<double>::infinity();
    //long double best_score = -std::numeric_limits<long double>::infinity();
    int best_j = -1;

    // Ensure j_start is valid
    int j_start = std::min(i - 1, ju);
    j_start = std::max(j_start, q - 1); // must leave enough points for left bins

    for (int j = j_start; j >= jl; --j) {
        // Each bin must have at least one point
        if ((i - j) < 1) continue;

        // Skip duplicate splits (optional)
        if (j > 0 && j != jl && fabs(data[j] - data[j - 1]) < 1e-12) continue;

        // double llh = loglikelihood(j + 1, i, index, entropy);
        
        double score = S[j][q - 1] + loglikelihood(j + 1, i, index, entropy);
        //long double score = static_cast<long double>(S[j][q - 1]) + loglikelihood(j + 1, i, index, entropy);

        // const double tol = 1e-10; // or 1e-14
        
        // cout << std::fixed << std::setprecision(20) << "score:\t" << score << "\t best_score:\t" << best_score <<
        //     "\t S:\t" << S[j][q - 1] << "\t llh:\t" << llh << endl;
        if (score >= best_score) {
            best_score = score;
            best_j = j;
        }
        // if (score > best_score || (fabs(score - best_score) < 1e-16 && j > best_j)) {
        //     best_score = score;
        //     best_j = j;
        // }
    }

    // Fallback if no valid j found
    if (best_j < 0) {
        best_j = jl;
        //best_score = S[best_j][q - 1] + loglikelihood(best_j + 1, i, index, entropy);
        best_score = (S[best_j][q - 1]) + loglikelihood(best_j + 1, i, index, entropy);
    }

    //S[i][q] = best_score;
    S[i][q] = static_cast<double>(best_score);
    P[i][q] = best_j;

    if (best_j >= 0) {
        cal_score_binomial_recur(q, il, jl, i - 1, best_j, S, P, index, data, entropy);
        cal_score_binomial_recur(q, i + 1, best_j, iu, ju, S, P, index, data, entropy);
    }
}

double loglikelihood(int start, int end, const vector<vector<int>>& index, bool entropy) {
    unsigned int total_points = 0;
    double score_region = 0.0;
    //long double score_region = 0.0L;
    size_t num_clusters = index.size();

    for (size_t cluster = 0; cluster < num_clusters; ++cluster) {
        unsigned int count = (start <= 0) ? index[cluster][end] : index[cluster][end] - index[cluster][start - 1];
        if (count == 0) continue;

        total_points += count;
        //score_region += count * std::log(static_cast<double>(count));
        score_region += count * std::log(count);
        // cout << std::fixed << std::setprecision(16) << "count:\t" << count <<
        //     "\t log(count):\t" << std::log(count) << "\t score_region:\t" << score_region << endl;
    }

    if (total_points == 0) return 0.0;

    double log_total = std::log(static_cast<double>(total_points));
    //long double log_total = std::log(static_cast<long double>(total_points));
    if (entropy) {
        score_region -= log_total;
    } else {
        score_region -= total_points * log_total;
        //score_region -= static_cast<long double>(total_points) * log_total;
    }
    // cout << std::fixed << std::setprecision(16) << "total_points:\t" << total_points <<
    //     "\t log_total:\t" << log_total << "\t score_region:\t" << score_region << endl;
    return score_region;
}

/*void cal_score_binomial_recur_main(const vector<double>& data, const vector<int>& labels_compressed,
                                   const vector<int>& labels_weights,
                                   vector<vector<double>>& S, vector<vector<int>>& P,
                                   const vector<vector<int>>& index,
                                   unsigned int Q_max, unsigned int K, int n_blocks, bool entropy) {
    // Initialize base case for q = 0 (first bin only)
    vector<int> init_weights(K, 0);
    vector<int> total_weights(n_blocks, 0);

    // maybe can initialize index here to avoid unnecessary traversal (resize)
    for (int i = 0; i < n_blocks; ++i) {
        init_weights[labels_compressed[i]] += labels_weights[i];
        total_weights[i] = (i == 0) ? labels_weights[i] : total_weights[i - 1] + labels_weights[i];


        double score = 0;
        for (size_t x = 0; x < K; ++x) {
            if (init_weights[x] > 0) {
                score += init_weights[x] * log(static_cast<double>(init_weights[x]) / total_weights[i]);
            }
        }
        S[i][0] = score;
        P[i][0] = i;
    }
    // Use divide-and-conquer DP optimization for q > 0
    for (size_t q = 1; q < Q_max; ++q) {
        cal_score_binomial_recur(q, q - 1, q - 1, n_blocks, n_blocks - 1, S, P, index, data, entropy);
    }
}

void cal_score_binomial_recur(int q, int il, int jl, int iu, int ju, vector<vector<double>>& S, vector<vector<int>>& P,
                              const vector<vector<int>>& index, const vector<double>& data, bool entropy) {
    if (il >= iu) return;

    int i = (il + iu) / 2;
    double best_score = -numeric_limits<double>::infinity();
    int best_j = jl;

    // Find best split point j for current i
    for (int j = min(i - 1, ju); j >= jl; --j) {
        if (j > 0 && data[j] == data[j - 1]) continue;

        double score = S[j][q - 1] + loglikelihood(j, i, index, entropy);
        if (score > best_score) {
            best_score = score;
            best_j = j;
        }
    }

    // Update DP table
    S[i][q] = best_score;
    P[i][q] = best_j;

    // Recursively solve left and right subproblems
    cal_score_binomial_recur(q, il, jl, i, best_j, S, P, index, data, entropy);
    cal_score_binomial_recur(q, i + 1, best_j, iu, ju, S, P, index, data, entropy);
}

double loglikelihood(int i, int j, const vector<vector<int>>& index, bool entropy) {
    unsigned int total_points = 0;
    double score_region = 0.0;
    const size_t num_clusters = index.size();

    // Compute total points and log-likelihood contribution from each cluster
    for (size_t cluster = 0; cluster < num_clusters; ++cluster) {
        unsigned int count = (i == -1) ? index[cluster][j] : index[cluster][j] - index[cluster][i];
        if (count == 0) continue;

        total_points += count;
        score_region += count * std::log(static_cast<double>(count));
    }

    if (total_points == 0) return 0.0;

    double log_total = std::log(static_cast<double>(total_points));
    if (entropy) {
        score_region -= log_total;
    } else {
        score_region -= total_points * log_total;
    }

    return score_region;
}*/

template <typename T>
void print_out_table(const vector<vector<T>>& vec, int width, int precision) {
    cout << "     "; // For alignment of column headers
    // Print column headers
    for (size_t j = 0; j < vec[0].size(); ++j) {
        cout << setw(width) << j << " ";
    }
    cout << endl;
    for (size_t i = 0; i < vec.size(); ++i) {
        cout << i << ": ";
        for (const auto& elem : vec[i]) {
            cout << setw(width) << fixed << setprecision(precision) << elem << " ";
        }
        cout << endl;
    }
}
