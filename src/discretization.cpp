//
// Created by Jiandong Wang on 2022/11/12.
//

#include <R.h>
#include "discretization.h"

//static const vector<double> *spx;

mutex m2;

/*static bool compi(size_t i, size_t j) {
    return (*spx)[i] < (*spx)[j];
}*/


vector<vector<int>>
discretize_data(Cluster &cluster, vector<vector<double>> lines_all_dim, size_t n_thread, size_t step, bool skip) {
    /**
     * @param lines_all_dim the separate lines for all dimension
     * @return a two dimension table, with col as dimension and row is sample. The table is a collection of index of
     * each point for all dimension. A row in the table store the index(label) for that dimension(row).
     */
    size_t ndim = cluster.get_dims();
    size_t nsample = cluster.get_sample_size();
    vector<vector<double>> data = cluster.cal_order();
    auto order = cluster.get_order();
    vector<vector<int>> discr_data(nsample, vector<int>(ndim, 0));
    size_t incre = 0;

/*    for (size_t dim = 0; dim < ndim; ++dim) {
        // loop for dims
        vector<double> lines = lines_all_dim[dim];
        // data on that dimension need to be sorted

        if (lines.size() == 0 and skip) {
            // no cut will be performed
            incre += 1;
            continue;
        } else {

            //create an order indicator
            vector<size_t> order(nsample);
            for (size_t i = 0; i < order.size(); ++i) {
                order[i] = i;
            }
            vector<double> sorted_data(nsample);
            spx = &data[dim];
            //order the points and re-organize the label
            sort(order.begin(), order.end(), compi);
            for (size_t i = 0; i < order.size(); ++i) {
                sorted_data[i] = data[dim][order[i]];
            }

            lines.push_back(sorted_data[nsample - 1]);
            auto iter = sorted_data.begin();
            auto iter_discr = discr_data.begin();
            double line = lines[0];
            size_t label = 0;
            // lines should be in increasing order, may have to be checked
            // todo: test if sort is necessary
            for (int i = 0; i < nsample; ++i) {
                // loop for all data points
                if ((*iter) <= line) {
                    //on the left side of the line(smaller than equal) , not included
                    discr_data[order[i]][dim - incre] = label;
                } else {
                    // on the right side (greater than)
                    ++label;
                    line = lines[label];
                    discr_data[order[i]][dim - incre] = label;
                }
                ++iter;
            }
        }
    }*/

    vector<thread> threads(n_thread);
    size_t dim_l, dim_u;
    for (auto i = 0u; i < n_thread; ++i) {
        dim_l = step * i;
        dim_u = step * (i + 1);
        if (i == n_thread - 1) {
            dim_u = ndim;
        }
        //cout << "thread " << i << " created for discretization!" << endl;
        threads[i] = thread(discretization_data_thread, ref(discr_data), cref(data), cref(order),
                            cref(lines_all_dim), nsample, dim_u, dim_l, skip);
    }

    for (auto &thread: threads) {
        thread.join();
        //cout << "thread " << thread.get_id() << " joined for discretization!" << endl;
    }

    for (size_t i = 0u; i < nsample; ++i) {
        discr_data[i].resize(ndim - incre);
    }
    return discr_data;
}

void discretization_data_thread(vector<vector<int>> &discr_data, const vector<vector<double>> &data,
                                const vector<vector<double>> &orders, const vector<vector<double>> &lines_all_dim,
                                size_t nsample, size_t dim_u, size_t dim_l, bool skip) {
    for (size_t dim = dim_l; dim < dim_u; ++dim) {
        // loop for dims
        vector<double> lines = lines_all_dim[dim];
        // data on that dimension need to be sorted

        if (lines.size() == 0 and skip) {
            // no cut will be performed
            //incre += 1;
            continue;
        } else {
            //create an order indicator
            auto order = orders[dim];

            //my_compare mc(data, order);
            //mc.order.resize(nsample);
/*            for (size_t i = 0; i < nsample; ++i) {
                order[i] = i;
                mc.order[i] = i;
            }*/
            vector<double> sorted_data(nsample);

/*            mc.spx = &mc.data[dim];
            mc.dim = dim;

            m2.lock();
            spx = &data[dim];
            order the points and re-organize the label
            sort(order.begin(), order.end(), mc.compi());
            sort(order.begin(), order.end(), compi);
            m2.unlock();*/
            //mc.sort_on_dim(dim);

            for (size_t i = 0; i < order.size(); ++i) {
                sorted_data[i] = data[dim][order[i]];
            }

            lines.push_back(sorted_data[nsample - 1]);
            auto iter = sorted_data.begin();
            // auto iter_discr = discr_data.begin();
            double line = lines[0];
            size_t label = 0;
            sort(lines.begin(), lines.end());
            for (size_t i = 0; i < nsample; ++i) {
                // loop for all data points
                if ((*iter) <= line) {
                    //on the left side of the line(smaller than equal) , not included
                    discr_data[order[i]][dim] = label;
                } else {
                    // on the right side (greater than)
                    ++label;
                    line = lines[label];
                    discr_data[order[i]][dim] = label;
                }
                ++iter;
            }
        }
    }
}

void _test_discretize_data() {
    vector<vector<double>> data = {
            {1, 2, 4},
            {1, 2, 4},
            {2, 2, 4},
            {3, 2, 4},
            {8, 6, 4},
            {4, 3, 4},
            {5, 3, 6},
            {6, 3, 6},
            {6, 4, 6},
            {7, 5, 6}
    };
    vector<int> labels(data.size(), 1);
    vector<vector<double>> lines = {{1.5, 5, 8},
                                    {3,   6},
                                    {}};
    Cluster cluster(labels, data);
    vector<vector<int>> discr_data;
    discr_data = discretize_data(cluster, lines, 1, cluster.get_dims(), false);
    for (size_t i = 0; i < discr_data.size(); ++i) {
        for (auto d: discr_data[i]) {
            Rprintf("%d,", d); // cout << d << ",";
        }
        Rprintf("\n"); // cout << endl;
    }

}
