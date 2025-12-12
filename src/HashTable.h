//
// Created by Jiandong Wang on 2023/5/1.
//

#ifndef TELESCOPE_HASHTABLE_H
#define TELESCOPE_HASHTAB LE_H

#include <unordered_map>
#include <vector>
#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>

using namespace std;

struct KeyForTele {
    size_t index; // index to an original data point. Data is not copied here to save memory
    KeyForTele(size_t i) : index(i) {}
};

template<>
struct std::hash<KeyForTele> {
    // Custom hash function class
    const vector<vector<size_t> > &m_data;
    vector<size_t> m_hash_vars;

    hash(const vector<vector<size_t> > &data) : m_data(data) {
        auto n = data.size();
        auto d = data[0].size();
        auto m = floor(log2(1 + n));
        m = m < d ? m : d;
        m_hash_vars.resize(m);
        for (auto j = 0; j < m; j++) {
            m_hash_vars[j] = j;
        }
    }

    std::size_t operator()(KeyForTele const &key) const noexcept {
        size_t h = 0;
        for (auto &v: m_hash_vars) {
            h ^= std::hash<int>{}(m_data[key.index][v]);
        }
        return h;
    }
};

template<>
struct std::equal_to<KeyForTele> {
    // Custom equal_to class

    const vector<vector<size_t> > &m_data;

    equal_to(const vector<vector<size_t> > &data) : m_data(data) {}

    bool operator()(const KeyForTele &x, const KeyForTele &y) const {
        return m_data[x.index] == m_data[y.index];
    }
};

struct CountForTele {
    size_t tot;
    // No Label needed
    vector<size_t> by_index;
    size_t posi = 0;

    CountForTele() {}

//    Count(size_t c) {}

//    Count(size_t c, size_t nindexs)
//            : tot(c), by_index(nindexs, 0) {}

    CountForTele(size_t c, size_t n_points, size_t index)
            : tot(c), by_index(n_points, 0) {
        by_index[posi] = index + 1;
        posi++;
    }
};

struct HashTableForTele {
    const vector<vector<size_t> > &m_data;
    std::hash<KeyForTele> m_hasher;
    std::equal_to<KeyForTele> m_key_equal;
    // m_table must be declared after m_hasher and m_key_equal:
    std::unordered_map<KeyForTele, CountForTele> m_table;
//    vector<size_t> groups;
//    size_t max_group = 0;

      HashTableForTele(const vector<vector<size_t> > &data)
            : m_data(data), m_hasher(data),
              m_key_equal(data), m_table(0, m_hasher, m_key_equal) {
        auto n = data.size();

        for (size_t i = 0u; i < n; ++i) {
            KeyForTele key(i);
            auto p = m_table.find(key);
            if (p == m_table.end()) {
               CountForTele c(1, n, i);
                m_table.insert({key, c});
            } else {
                p->second.tot++;
                p->second.by_index[p->second.posi] = i + 1;
                p->second.posi++;
            }
        }
    }

/*    Count operator()(size_t i) noexcept {
        return m_table[Key(i)];
    }*/

    void print() {
        for (auto & p : m_table) {
            auto &v = m_data[(p.first).index];
            for (auto &x: v) cout << x << ",";
            cout << ": " << p.second.tot << ";\t";
            //auto &l = p->second.by_index;
            for (size_t x = 0; x < p.second.posi; ++x) {
                cout << p.second.by_index[x] << ",";
            }
            cout << ";\t" << p.second.posi;
            cout << endl;
        }
    }

    vector<size_t> gen_label() {
        vector<size_t> labels(m_data.size(), 0);
        int i = 0;
        for (auto & p : m_table) {
            for (size_t x = 0; x < p.second.posi; ++x) {
                labels[p.second.by_index[x]-1] = i;
            }
            ++i;
        }
        return labels;
    }

/*    std::vector<vector<size_t>> gen_table() {
        size_t nlabels = 1 + *std::max_element(m_labels.begin(), m_labels.end());
        std::vector<vector<size_t>> conti_table(m_table.size(), std::vector<size_t>(nlabels, 0));
        auto iter_conti = conti_table.begin();

        // cluster is row, first loop
        for (auto p = m_table.begin(); p != m_table.end(); ++p) {
            // *iter_conti = p->second.by_label;
            iter_conti++;
        }

        for (auto i: conti_table) {
            for (auto j: i) {
                cout << j << ",";
            }
            cout << endl;
        }
        return conti_table;
    }*/
};

#endif //TELESCOPE_HASHTABLE_H
