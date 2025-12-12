/**
 * @file Joint_Grid.cpp
 * @author Jiandong Wang
 * @date 2/15/20
 * @copyright NMSU Song lab
 */

#include <R.h>
#include "Joint_Grid.h"
#include <cfenv>
// #pragma STDC FENV_ACCESS ON

mutex m0, m1;

struct grid Find_Grid(
    Cluster& clusters, 
    const vector<int>& min_bin_limit,
    const vector<int>& max_bin_limit, 
    const string& method,
    double cut_off, bool entropy, 
    bool reduction, size_t n_thread
)
{
    std::fesetround(FE_TONEAREST);
    /**
     *  @param clusters a cluster object contains the points and class information of each cluster.
     *  @param min_bin_limit a vector of integers that stores the minimum number of bins for each dimension.
     *  @param max_bin_limit a vector of integers that stores the maximum number of bins for each dimension.
     *  @param method a string passing the discretization method user selected. Should be one of the following:
     *  "Density":
     *  "DP majority": ,
     *  "DP compressed majority":
     *  "DP exact likelihood":
     *  "DP approx likelihood":
     *  @param cut_off
     *  @param entropy
     *  @details A discretization function for numerical data, especially designed for high dimensional data.
     *
     */
    size_t n_dims = clusters.get_dims();
    size_t N = clusters.get_sample_size();
    size_t MAX_thread = thread::hardware_concurrency();
    if (MAX_thread == 0) {
        MAX_thread = 1; // fallback to single thread
    }

    // checking and prepare for multi-thread
    if (n_thread <= 0 || n_thread > MAX_thread || n_thread > n_dims)
    {
      // cout << 
      REprintf(
        "number of threads error, automatically adapted!\n");
      n_thread = min(MAX_thread, n_dims);
    }
    //cout << "number of thread: " << n_thread << endl;
    size_t step = n_dims / n_thread;
    //n_thread = n_dims / step + 1;

    clusters.cal_order();
    vector<vector<double>> lines_multi_dim(n_dims);
    vector<vector<double>> medians_multi_dim(n_dims);
    vector<size_t> num_lines(n_dims);
    vector<thread> threads(n_thread);
    vector<double> best_bic_all(n_dims);
    size_t dim_l, dim_u;
    for (auto i = 0u; i < n_thread; ++i)
    {
        dim_l = step * i;
        dim_u = step * (i + 1);
        if (i == n_thread - 1)
        {
            dim_u = n_dims;
        }
        //cout << "thread " << i << " created!" << endl;
        threads[i] = thread(Find_Grid_thread, ref(clusters), min_bin_limit, max_bin_limit, cref(method),
                            cut_off, entropy, ref(lines_multi_dim), ref(num_lines), ref(medians_multi_dim),
                            dim_u, dim_l, ref(best_bic_all));
        /*        vector<double> lines = Find_1D_Grid(clusters, i, min_bin_limit[i], max_bin_limit[i], method,
                                                    cut_off, entropy);
                lines_multi_dim[i] = lines;
                num_lines[i] = lines.size();*/
    }

    for (auto& thread : threads)
    {
        thread.join();
        //cout << "thread " << thread.get_id() << " joined!" << endl;
    }

    grid G;
    G.lines = lines_multi_dim;
    G.num_lines = num_lines;
    G.dimensional_bic = best_bic_all;

    //discretization
    vector<vector<int>> discr_data = discretize_data(clusters, lines_multi_dim, n_thread, step, false);

    //hash the table into hashtable

    vector<size_t> labels_t(N, 0);
    auto iter_lab = labels_t.begin();
    for (int label : clusters.get_labels())
    {
        *iter_lab = abs(label);
        ++iter_lab;
    }

    MyHash hct(discr_data, labels_t);

    //ct.print();
    //Conti-table: row is the cluster(predict), col is the class(truth),

    auto conti_table = hct.gen_table();

    double BIC_table = BIC_for_table(conti_table, N, num_lines);
    double purity_original = purity(conti_table, N);

    G.conti_table = conti_table;
    G.original_BIC = BIC_table;
    G.original_purity = purity_original;
    G.discr_data = discr_data;

    // calculate upsilon score

    auto upsilon = upsilon_stat(conti_table);
    G.upsilon = upsilon;

    if (reduction)
    {
        /*
                vector<bool> dim_removed(n_dims, false);
                vector<vector<double>> line_after_removed = lines_multi_dim;
                vector<vector<size_t>> best_conti_table = conti_table;
                double BIC_after_removed = 0;

                dim_removed = dim_reduct(clusters, line_after_removed, best_conti_table, BIC_table, BIC_after_removed, n_thread,
                                         step);

                G.lines_after_removed = line_after_removed;
                G.line_removed = dim_removed;
                G.conti_table_after = best_conti_table;
                G.removed_BIC = BIC_after_removed;
                //discretize_data(clusters, lines_multi_dim);
                G.discr_data_after = discretize_data(clusters, line_after_removed, n_thread, step);
                //todo: calculate ARI and Purity

                double purity_after = purity(best_conti_table, N);
                G.removed_purity = purity_after;*/
    }
    else
    {
        G.removed_BIC = 0;
        G.removed_purity = 0;
        G.line_removed = vector<bool>(n_dims, false);
    }

    vector<vector<double>>  median_distances=vector<vector<double>>(hct.m_table.size(), vector<double>(N, 0));
    vector<vector<double>>  mean_distances=vector<vector<double>>(hct.m_table.size(), vector<double>(N, 0));

    auto all_points = clusters.get_all_points();

    // auto iter_dist = distances.begin();
    // for (auto p = hct.m_table.begin(); p != hct.m_table.end(); ++p)
    // {
    //     auto data_sub = getSubsetByIndices(all_points, p->second.members);
    //     vector<double> median = computeMedians(data_sub);
    //     auto temp_dist = computeDistances(data_sub, median);
    //     *iter_dist = temp_dist;
    //     ++iter_dist;
    // }

    median_distances = evaluation_Distances(all_points, hct, "median");
    mean_distances = evaluation_Distances(all_points, hct, "mean");

    G.median_distance = median_distances;
    G.mean_distance = mean_distances;

    clusters.set_grids(G);

    return G;
}

///@param clusters
///@param dim_l lower bound of the dimension to be processed, included
///@param dim_u upper bound of the dimension to be processed, not included
void Find_Grid_thread(Cluster& clusters, const vector<int>& min_bin_limit, const vector<int>& max_bin_limit, const string& method,
                      double cut_off, bool entropy, vector<vector<double>>& lines_multi_dim, vector<size_t>& num_lines,
                      vector<vector<double>>& medians_multi_dim, size_t dim_u, size_t dim_l, vector<double> & best_bic_all)
{
    for (auto i = dim_l; i < dim_u; ++i) {
        double selected_bic = numeric_limits<double>::infinity();
        vector<double> lines;
        if(method == "DP approx likelihood two way"){
          double best_bic_inc = numeric_limits<double>::infinity();
          double best_bic_dec = numeric_limits<double>::infinity();

          vector<double> lines_inc = Find_1D_Grid(clusters, i, min_bin_limit[i], max_bin_limit[i], medians_multi_dim[i],
                                                  method, cut_off, entropy, best_bic_inc, true);

          vector<double> lines_dec = Find_1D_Grid(clusters, i, min_bin_limit[i], max_bin_limit[i], medians_multi_dim[i],
                                                  method, cut_off, entropy, best_bic_dec, false);
            // Pick direction with *lower* BIC (better)
            bool choose_inc = (best_bic_inc <= best_bic_dec);

            if (choose_inc) {
                //cout << "Select Increasing" << endl;
                lines = std::move(lines_inc);
                selected_bic = best_bic_inc;
            } else {
                //cout << "Select Decreasing" << endl;
                std::sort(lines_dec.begin(), lines_dec.end());
                lines = std::move(lines_dec);
                selected_bic = best_bic_dec;
            }
       }else {
           lines = Find_1D_Grid(clusters, i, min_bin_limit[i], max_bin_limit[i], medians_multi_dim[i],
                                                           method, cut_off, entropy, selected_bic, true);
        }

        {
            std::lock_guard<std::mutex> lock(m0);
            lines_multi_dim[i] = std::move(lines);
            num_lines[i] = lines_multi_dim[i].size();
            best_bic_all[i] = selected_bic;
        }
    }
}


vector<double>
Find_1D_Grid(Cluster& clusters, int dim_input, int min_bin, int max_bin, vector<double>& medians, const string& method,
             double cut_off, bool entropy, double & best_bic, const bool& increasing)
    {
    if (method == "Density")
    {
        if (not clusters.split_by_cluster())
        {
            clusters.split_data();
        }
        if (not clusters.has_median())
        {
            clusters.cal_medians();
        }

        return (Cut_by_density_1D(clusters, dim_input, min_bin));
    }
    else if (method == "DP majority" or method == "DP compressed majority" or method == "DP exact likelihood" or
        method == "DP approx likelihood" or method == "DP approx likelihood two way")
    {
        /**
         * @note methods name changed:\n
         * dp compressed likelihood -> dp exact likelihood\n
         * dp compressed likelihood nlogn -> dp approx likelihood
         */
        // order the point to get the list of label
        int K = clusters.get_K();
        int N = clusters.get_sample_size();

        vector<double> points_all = clusters.get_all_points_by_dim(dim_input);
        vector<int> labels_all = clusters.get_labels();
        vector<double> points_all_sorted(N);
        vector<int> labels_all_sorted(N);
        /*        //cat all clusters together
                for (int k = 0; k < K; ++k) {
                    vector<int> labels_temp(points_dim[k].size(), k);
                    points_all.insert(points_all.end(), points_dim[k].begin(), points_dim[k].end());
                    labels_all.insert(labels_all.end(), labels_temp.begin(), labels_temp.end());
                }*/

        if (test_sorted(points_all, increasing))
        {
            points_all_sorted = points_all;
            labels_all_sorted = labels_all;
        }
        else
        {
            clusters.cal_order();
            auto order = clusters.get_order()[dim_input];

            // If backward sorting is needed, reverse the order vector
            if (!increasing) {
                std::reverse(order.begin(), order.end());
            }

            for (size_t i = 0; i < order.size(); ++i)
            {
                points_all_sorted[i] = points_all[order[i]];
                labels_all_sorted[i] = labels_all[order[i]];
            }
        }

        vector<int> loca;
        if (method == "DP majority")
        {
            loca = Cutting_Cluster(labels_all_sorted, K, max_bin, min_bin, cut_off);
        }
        else if (method == "DP compressed majority" or method == "DP exact likelihood" or
            method == "DP approx likelihood" or method == "DP approx likelihood two way")
        {
            loca = Cutting_Cluster_dp_compressed(labels_all_sorted, points_all_sorted, K, max_bin,
                                                 min_bin, best_bic, method, entropy);
        }


        /*vector<double> lines;(loca.size() - 1);
        for (size_t i = 0; i + 1 < loca.size(); ++i)
        {
            lines[i] = (points_all_sorted[loca[i]] + points_all_sorted[loca[i] + 1]) / 2.0;
        }*/

        vector<double> lines;
        if(loca.size()>0)
        {
            lines.resize(loca.size() - 1);
            for (size_t i = 0; i + 1 < loca.size(); ++i)
            {
                lines[i] = (points_all_sorted[loca[i]] + points_all_sorted[loca[i] + 1]) / 2.0;
            }
        }
        else
        {
            lines.push_back(points_all_sorted.back());
        }

        /*vector<double> median(loca.size());
        for (size_t i = 0; i < loca.size(); ++i)
        {
            int range = (i == 0) ? loca[i] + 1 : loca[i] - loca[i - 1];
            int midIndex = (i == 0) ? loca[i] / 2 : ((loca[i] + loca[i - 1] + 1) / 2);

            if (range % 2 == 0)
            {
                median[i] = (points_all_sorted[midIndex] + points_all_sorted[midIndex + 1]) / 2;
            }
            else
            {
                median[i] = points_all_sorted[midIndex];
            }
        }
        medians = median;*/
        return lines;
    }
    else
    {
        assert("method does not exist!");
    }
    // return {0};

    return {};
}

bool test_sorted(const vector<double>& x, bool increasing)
{
    if (x.size() < 2) return true;

    if (increasing) {
        for (size_t i = 0; i < x.size() - 1; ++i)
            if (x[i] > x[i + 1]) return false;
    } else {
        for (size_t i = 0; i < x.size() - 1; ++i)
            if (x[i] < x[i + 1]) return false;
    }

    return true;
}


double euclidean_distance(const vector<double>& point1, const vector<double>& point2)
{
    if (point1.size() != point2.size()) { throw invalid_argument("Points must have the same dimension"); }
    double sum = 0.0;
    for (size_t i = 0; i < point1.size(); ++i) { sum += pow(point1[i] - point2[i], 2); }
    return sqrt(sum);
}
