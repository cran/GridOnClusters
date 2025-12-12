//
// Created by Jiandong Wang on 2022/5/8.
//

#include "evaluation.h"

using namespace std;

// Choose the optimal number of levels between Kmin and Kmax for entropy approach
double evaluation_entropy(double N, double Q, double K, vector<vector<double>> S, vector<vector<double>> S_entropy) {
    ///\b clean version of the select_level function by removing all unnecessary part\n
    ///\b N is the number of samples (points)\n
    ///\b Q is the number of bins\n
    ///\b K is the number of groups
    ///clusters (channels)\n using C in 1-D case\n
    ///as we use the pre-calculated data (likelihood and entropy) from dynamic programming step, no calculation needed.

    if (N < 2) {
        return Q;
    }

    double log_likelihood = S[S.size() - 1][Q - 1];
    double bic = -2 * log_likelihood + (Q - 1 + (K - 1) * Q) * log((double) N);
    double bic_h = bic - N * S_entropy[S.size() - 1][Q - 1];
    return bic_h;
}

// Choose the optimal number of levels between Kmin and Kmax for likelihood method
double evaluation_likelihood(double N, double Q, double K, vector<vector<double>> S) {
    ///\b clean version of the select_level function by removing all unnecessary part\n
    ///\b N is the number of samples (points)\n
    ///\b Q is the number of bins\n
    ///\b K is the number of groups
    ///clusters (channels)\n using C in 1-D case\n
    /// in this version, as we use the pre-calculated data (likelihood) from dynamic programming step, no calculation needed.

    if (N < 2) {
        return Q;
    }

    double log_likelihood = S[S.size() - 1][Q - 1];
    double bic = -2 * log_likelihood + (Q - 1 + (K - 1) * Q) * log((double) N);
    return bic;
}

// Choose the optimal number of levels between Kmin and Kmax
double
evaluation_counts(const vector<double> &points, const vector<int> &labels, vector<vector<int>> index, vector<int> posi,
                  vector<vector<int>> &table) {
    ///N is the number of samples (points)\n
    ///Q is the number of bins\n
    ///K is the number of clusters (channels)\n
    /// in this version, as we use the pre-calculated data (likelihood) from dynamic programming step, no calculation needed.
    const size_t N = points.size();
    const size_t Q = posi.size();
    const size_t K = table.size();
    vector<vector<double>> channels = split_into_channels(points, labels, K, N);
    if (N < 2) {
        return (double) Q;
    }
    vector<int> size(Q);

    /// calculate the size
    size[0] = posi[0] + 1;
    for (size_t i = 1; i < Q; ++i) {
        size[i] = posi[i] - posi[i - 1];
    }

    vector<double> lambda(Q);
    vector<double> mu(Q);
    vector<double> sigma2(Q);
    vector<double> coeff(Q);

    size_t indexLeft = 0;
    size_t indexRight;

    /// Estimate GMM parameters for pool
    for (size_t q = 0; q < Q; ++q) {
        lambda[q] = size[q] / (double) N;

        indexRight = indexLeft + size[q] - 1;

        shifted_data_variance(points, indexLeft, indexRight, mu[q], sigma2[q]);

        if (sigma2[q] == 0 || size[q] == 1) {

            double dmin;

            if (indexLeft > 0 && indexRight < N - 1) {
                dmin = min(points[indexLeft] - points[indexLeft - 1], points[indexRight + 1] - points[indexRight]);
            } else if (indexLeft > 0) {
                dmin = points[indexLeft] - points[indexLeft - 1];
            } else {
                dmin = points[indexRight + 1] - points[indexRight];
            }

            if (sigma2[q] == 0) sigma2[q] = dmin * dmin / 4.0 / 9.0;
            if (size[q] == 1) sigma2[q] = dmin * dmin;
        }


/*
        if(sigma2[q] == 0) sigma2[q] = variance_min;
        if(size[q] == 1) sigma2[q] = variance_max;
*/

        coeff[q] = lambda[q] / sqrt(2.0 * M_PI * sigma2[q]);

        indexLeft = indexRight + 1;
    }

    double loglikelihood = 0;
/*
    vector<int> bin_sums(table[0].size(), 0);
    vector<int> bin_max(table[0].size(), 0);
*/

/*    /// calculate the sum and the max of bins
   for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            bin_sums[i] += table[j][i];
            if (bin_max[i] < table[j][i]) {
                bin_max[i] = table[j][i];
            }
        }
        if (bin_sums[i] == 0 or (bin_sums[i] == bin_max[i])) {
            continue;
        }
        loglikelihood *= (pow(((double) bin_max[i] / bin_sums[i]), bin_max[i]) *
                          pow(((double) (bin_sums[i] - bin_max[i]) / bin_sums[i]), (bin_sums[i] - bin_max[i])));
        loglikelihood += (bin_max[i] * log((double) bin_max[i] / bin_sums[i])) +
                         (bin_sums[i] - bin_max[i]) * log((double) (bin_sums[i] - bin_max[i]) / bin_sums[i]);
    }

    loglikelihood = log(loglikelihood);
*/
    for (size_t i = 0; i < N; ++i) {
        double L = 0;
        for (size_t k = 0; k < Q; ++k) {
            L += coeff[k] * exp(-(points[i] - mu[k]) * (points[i] - mu[k]) / (2.0 * sigma2[k]));
        }
        loglikelihood += log(L);
    }

    vector<double> loglikelihood_k(Q, 0);
    /// Estimate GMM parameters for channels
    for (size_t c = 0; c < K; ++c) {
        vector<double> lambda_c(Q);
        vector<double> mu_c(Q);
        vector<double> sigma2_c(Q);
        vector<double> coeff_c(Q);
        //vector <vector<double>> coeff_c(K, vector<double>(K));

        ///calculare the size for each channel
        vector<int> size_c(Q);
        size_c[0] = index[c][0];
        for (size_t i = 1; i < Q; ++i) {
            size_c[i] = index[c][i] - index[c][i - 1];
        }

        //size_t indexLeft = 0;
        //size_t indexRight = 0;

        for (size_t q = 0; q < Q; ++q) { // Estimate GMM parameters for each q
            lambda_c[q] = size_c[q] / (double) channels[c].size();

            indexRight = indexLeft + size_c[q] - 1;

            //empty bin (lambda is 0), skip, mean = 0, sigma = 1

            shifted_data_variance(channels[c], indexLeft, indexRight, mu_c[q], sigma2_c[q]);

            if (sigma2_c[q] == 0 || size_c[q] == 1) {

                double dmin;

                if (indexLeft > 0 && indexRight + 1 < channels[c].size()) {
                    dmin = min(channels[c][indexLeft] - channels[c][indexLeft - 1],
                               channels[c][indexRight + 1] - channels[c][indexRight]);
                } else if (indexLeft > 0) {
                    dmin = channels[c][indexLeft] - channels[c][indexLeft - 1];
                } else {
                    dmin = channels[c][indexRight + 1] - channels[c][indexRight];
                }

                if (sigma2_c[q] == 0) sigma2_c[q] = dmin * dmin / 4.0 / 9.0;
                if (size_c[q] == 1) sigma2_c[q] = dmin * dmin;
            }

            coeff_c[q] = lambda_c[q] / sqrt(2.0 * M_PI * sigma2_c[q]);

            indexLeft = indexRight + 1;
        }

        for (size_t i = 0; i < channels[c].size(); ++i) {
            double L = 0;
            for (size_t k = 0; k < Q; ++k) {
                L += coeff_c[k] * exp(-(channels[c][i] - mu_c[k]) * (channels[c][i] - mu_c[k]) / (2.0 * sigma2_c[k]));
            }
            loglikelihood_k[c] += log(L);
        }
    }

    vector<vector<int>> size_channels(K, vector<int>(Q, 0));
    vector<vector<double>> loglikelihood_c_k(K, vector<double>(Q, 0));
    /// without Estimate GMM parameters for channels
    for (size_t c = 0; c < K; ++c) {
        ///calculare the size for each channel
        vector<int> size_c(Q);
        size_c[0] = index[c][0];
        for (size_t i = 1; i < Q; ++i) {
            size_c[i] = index[c][i] - index[c][i - 1];
        }
        size_channels[c] = size_c;

        // int indexLeft = 0;
        // int indexRight = 0;

        for (size_t i = 0; i < channels[c].size(); ++i) {
            double L = 0;
            for (size_t k = 0; k < Q; ++k) {
                //L += coeff[k] * exp(-(channels[c][i] - mu[k]) * (channels[c][i] - mu[k]) / (2.0 * sigma2[k]));
                L += exp(-(channels[c][i] - mu[k]) * (channels[c][i] - mu[k]) / (2.0 * sigma2[k]));
            }
            loglikelihood_k[c] += log(L);
        }

        for (size_t i = 0; i < channels[c].size(); ++i) {
            int k = 0;
            for (size_t j = 0; j < index[c].size(); ++j) {
                if (i < (size_t) index[c][j]) {
                    k = j;
                    break;
                }
            }

            /// likelihood for each bin
            double L = 0;
            for (size_t q = 0; q < Q; ++q) {
                L += coeff[q] * exp(-(channels[c][i] - mu[q]) * (channels[c][i] - mu[q]) / (2.0 * sigma2[q]));
                //L += (1 / sqrt(2.0 * M_PI * sigma2[q])) *
                //     exp(-(channels[c][i] - mu[q]) * (channels[c][i] - mu[q]) / (2.0 * sigma2[q]));
            }
            L = (1 / sqrt(2.0 * M_PI * sigma2[k])) *
                exp(-(channels[c][i] - mu[k]) * (channels[c][i] - mu[k]) / (2.0 * sigma2[k]));
            loglikelihood_k[k] += log(L);

            /// likelihood for each channel each bin
            //double L = coeff[k] * exp(-(channels[c][i] - mu[k]) * (channels[c][i] - mu[k]) / (2.0 * sigma2[k]));
            //double L_c =(1 / sqrt(2.0 * M_PI * sigma2[k])) *  exp(-(channels[c][i] - mu[k]) * (channels[c][i] - mu[k]) / (2.0 * sigma2[k]));
            loglikelihood_c_k[c][k] += log(L);
        }
    }

    /* {
      /// sum of likelihood for each channels
      // double log_L = 0;
      /// calculate dLL
      for (size_t q = 0; q < Q; ++q) {
        // double log_L_c = 0;
        /// calculate dll(q)
        for (size_t m = 0; m < K; ++m) {
          // new equation
          // dLL(q) = 1/K sum_{j=1}^K  sum_{m=1}^K  | [n(m,q)*log (n(m,q)/n(m)) + L_m(q)] - [ n(j,q) * log (n(j,q)/n(j)) + L_j(q) ]  |
          for (size_t j = 0; j < K; ++j) {
            double pram_1, pram_2;
            if (size_channels[m][q] == 0) {
              pram_1 = 0;
            } else {
              pram_1 = size_channels[m][q] * log((double) size_channels[m][q] / channels[m].size());
            }
            if (size_channels[j][q] == 0) {
              pram_2 = 0;
            } else {
              pram_2 = size_channels[j][q] * log((double) size_channels[j][q] / channels[j].size());
            }
            // log_L_c += abs((pram_1 + loglikelihood_c_k[m][q]) - (pram_2 + loglikelihood_c_k[j][q]));
          }
          // log_L_c += abs(loglikelihood_k[q] - K * loglikelihood_c_k[c][q]);
        }
        
        // log_L += 1.0 / K * log_L_c;
      }
      
    } */


    /// Compute the Bayesian information criterion
    // double multi = 1;
    // for (size_t i = 0; i < Q; ++i) {
    //     multi *= size[i];
    // }
    //loglikelihood = S[S.size() - 1][K - 1];
    double bic = -2 * loglikelihood + (Q - 1 + (K - 1) * Q) * log((double) N);
    /*
    //double bic = -2 * loglikelihood + (2 * K - 1) * log((double) N);
    //double bic = - 2 * loglikelihood + K * log((double) N);
    //double bic = loglikelihood - log_L + (3 * K - 1 + (K * (K - 1))) * log((double) N);  //(K*3-1)
    //double bic = -2 * bic_c + 2 * loglikelihood - (3 * K - 1) * log((double) N);  //(K*3-1)
*/
    return bic;
}


void shifted_data_variance(const vector<double> &x, int left, int right, double &mean,
                           double &variance) {
    /// this function only will be used when BIC is calculated by absolute counts
    double sum = 0.0;
    double sumsq = 0.0;

    mean = 0.0;
    variance = 0.0;

    int n = right - left + 1;

    if (right >= left) {

        double median = x[(left + right) / 2];

        for (int i = left; i <= right; ++i) {
            sum += x[i] - median;
            sumsq += (x[i] - median) * (x[i] - median);
        }
        mean = sum / n + median;

        if (n > 1) {
            variance = (sumsq - sum * sum / n) / (n - 1);
        }
    }

}

vector<vector<double>>
split_into_channels(const vector<double> &points, const vector<int> &labels, unsigned int K, unsigned int N) {
    vector<vector<double>> cluster_points(K, vector<double>(N, 0));
    vector<size_t> clusters_size(K, 0);
    auto label_iter = labels.begin();
    auto data_iter = points.begin();
    int label;
    //loop for clusters
    for (; label_iter != labels.end() and data_iter != points.end(); ++label_iter, ++data_iter) {
        label = *label_iter;
        cluster_points[label][clusters_size[label]] = (*data_iter);
        clusters_size[label]++;
    }
    for (size_t i = 0; i < K; ++i) {
        (cluster_points[i]).resize(clusters_size[i]);
    }
    return cluster_points;
}
