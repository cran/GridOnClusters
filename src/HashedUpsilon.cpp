// HashedUpsilon.cpp
//
//   Perform the Upsilon test on a hashed contingency table
//
// Created: April 2, 2023. MS

#include <R.h>
#include "HashedUpsilon.h"

void upsilon_test_by_hashing(
    const vector<vector<int> > & data,
    double & upsilon, double & df, double & effect)
{
  upsilon = 0;
  df = 0;
  effect = 0;

  auto n = data.size();

  if(n > 0) {

    HCountTable hct(data);

    auto D = data[0].size();

    // Obtain marginal counts for all dimensions
    vector<unordered_map<int, int> > margin_counts(D);

    for(auto p = hct.m_table.begin(); p != hct.m_table.end(); ++p) {
      auto & l = data[(p -> first).m_i];
      auto & count = (p -> second).tot;
      for(auto d=0u; d<D; ++d) {
        margin_counts[d][l[d]] += count;
      }
    }

    // Calculate upsilon 0, assumming all entries are zero
    upsilon = n;
    for(auto d=0u; d<D; ++d) {
      auto & nd = margin_counts[d];
      double stat = 0.0;
      for(auto p=nd.begin(); p != nd.end(); ++p) {
         double Pd = p -> second / (double) n; // marginal probability
         stat += Pd * Pd;
      }
      upsilon *= nd.size() * stat;
    }

    // Calculate test statistic for each non-empty entry
    //   in the hashed contingency table
    double stat = 0;
    for(auto p = hct.m_table.begin(); p != hct.m_table.end(); ++p) {
      auto & l = data[(p -> first).m_i];
      auto & count = (p -> second).tot;
      double expected = n;
      for(auto d=0u; d<D; ++d) {
        expected *= margin_counts[d][l[d]] / (double) n;
      }
      stat += count * (count - 2.0 * expected);
    }

    df = 0;
    double num_entries = 1.0; // can overflow on very high dimensions

    for(auto d=0u; d<D; ++d) {
      auto & nd = margin_counts[d];
      num_entries *= nd.size();
      df += nd.size() - 1;
    }

    upsilon += stat / n * num_entries;

    df = num_entries - (df + 1);

    effect = sqrt(upsilon /
      (n * num_entries / pow(2.0, D) * (pow(2.0, D-1)-1)));

  }
  return;
}

// @title Test Upsilon by Hashing
// NO @export
// [[Rcpp::export]]
void test_upsilon_by_hashing()
{
  vector<vector<int> > data = {
    {1, 3},
    {1, 3},
    {1, 7},
    {1, 3},
    {1, 7},
    {1, 7},
  };

  double stat, df, effect;

  double stat_truth = 0;
  upsilon_test_by_hashing(data, stat, df, effect);
  if(stat != stat_truth) {
    REprintf(
      "ERROR in %s: %s line %d", 
      __FILE__, __func__, __LINE__
    );
    
    // std::cerr << "ERROR in " << __FILE__ << __func__ <<
    //  " line " << __LINE__ << "!" << endl;
  }

  data = {
    {2, 3},
    {2, 3},
    {1, 7},
    {2, 3},
    {1, 7},
    {1, 7},
  };

  /*
  pval_truth = 0.01430587844;
  pval = upsilon_test_by_hashing(data);
  if(abs(pval - pval_truth) > 0.00000000001) {
    cerr << "ERROR in test_upsilon_by_hashing(): 107!" << endl;
  }
  */

  stat_truth = 6;
  upsilon_test_by_hashing(data, stat, df, effect);
  if(stat != stat_truth) {
    REprintf(
      "ERROR in %s: %s line %d", 
      __FILE__, __func__, __LINE__
    );
    
    // cerr << "ERROR in " << __FILE__ << __func__ <<
    //  " line " << __LINE__ << "!" << endl;
  }

  data = {
    {1, 3, 2, 4, 5},
    {0, 3, 0, 4, 0},
    {1, 7, 2, 0, 6},
    {0, 3, 0, 4, 2},
    {1, 7, 2, 0, 7},
    {0, 0, 2, 4, 6},
    {0, 0, 0, 0, 0}
  };

  data = {
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0},
      {1, 1, 1, 1},
      {1, 1, 1, 1},
      {1, 1, 1, 1},
      {1, 1, 1, 1}
    };

  /*
  pval_truth = 0.01430587844;
  pval = upsilon_test_by_hashing(data);
  cout << "P = " << pval << endl;
  if(abs(pval - pval_truth) > 0.00000000001) {
    cerr << "ERROR in test_upsilon_by_hashing(): 123!" << endl;
  }
  */

  upsilon_test_by_hashing(data, stat, df, effect);
  
  Rprintf("Upsilon = %f, df = %f, effect = %f\n",
          stat, df, effect);
  
  // std::cout << "Upsilon = " << stat << ", df = " << df <<
  //  ", effect = " << effect << std::endl;

/*
  stat_truth = 10;

  if(stat != stat_truth) {
    cerr << "ERROR in test_upsilon_by_hashing(): 107!" << endl;
  }
*/

}
