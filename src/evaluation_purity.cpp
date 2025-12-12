//
// Created by Jiandong Wang on 2022/12/5.
//

#include "evaluation.h"
#include <iostream>

using namespace std;

double purity(vector<vector<size_t>> table, size_t N) {
    size_t sum = 0;
    for (auto i: table) {
        sum += *max_element(i.begin(), i.end());
    }
    return double(sum) / N;
}
