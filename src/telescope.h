//
// Created by Jiandong Wang on 2023/4/29.
//

#ifndef TELESCOPE_TELESCOPE_H
#define TELESCOPE_TELESCOPE_H

#include <vector>
#include <string>
#include <cmath>
#include "HashTable.h"

std::vector<size_t> tele(const std::vector<std::vector<double>> &data, const size_t n_lines);

void split(const std::vector<std::vector<double>> &data, const size_t n_lines, std::vector<std::vector<double>> &lines,
           std::vector<std::vector<size_t>> &tags_num_by_dim, std::vector<std::vector<size_t>> &tags_num);

std::vector<std::vector<size_t>> split(const std::vector<std::vector<double>> &data, const size_t n_lines);

std::vector<std::vector<double>> transform_by_dim(const std::vector<std::vector<double>> &data);

//std::vector<size_t> gen_label(const HCountTable hct, size_t n);

std::vector<std::vector<size_t>> sort_data(const std::vector<std::vector<double>> &data);
#endif //TELESCOPE_TELESCOPE_H
