//
//  HCountTable.hpp
//  HashedContingencyTable
//
//  Created by Joe M Song on 7/25/22.
//

#ifndef HCountTable_hpp
#define HCountTable_hpp
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <functional>
#include <algorithm>

using namespace std;

struct Key {
    size_t m_i; // index to an original data point. Data is not copied here to save memory

    Key(size_t i) : m_i(i) {}
};

template<>
struct std::hash<Key> {
    // Custom hash function class

    const vector<vector<int> > &m_data;

    vector<size_t> m_hash_vars;

    hash(const vector<vector<int> > &data) : m_data(data) {
        auto n = data.size();
        auto d = data[0].size();

        auto m = floor(log2(1 + n));

        m = m < d ? m : d;

        m_hash_vars.resize(m);

        for (auto j = 0; j < m; j++) {
            m_hash_vars[j] = j;
        }
    }

    std::size_t operator()(Key const &key) const noexcept {
        size_t h = 0;
        for (auto &v: m_hash_vars) {
            h ^= std::hash<int>{}(m_data[key.m_i][v]);
        }
        return h;
    }
};

template<>
struct std::equal_to<Key> {
    // Custom equal_to class

    const vector<vector<int> > &m_data;

    equal_to(const vector<vector<int> > &data) : m_data(data) {}

    bool operator()(const Key &x, const Key &y) const {
        return m_data[x.m_i] == m_data[y.m_i];
    }
};

struct Count {
    size_t tot;
    vector<size_t> by_label;
    vector<int> members;

    Count() : tot(0), by_label(0) {}

    Count(size_t c) : tot(c), by_label(0) {}

    Count(size_t c, size_t nlabels)
            : tot(c), by_label(nlabels, 0) {}

    void accrue(const Count &c) {
        tot += c.tot;
        if (by_label.size() == c.by_label.size()) {
            for (auto j = 0u; j < by_label.size(); ++j) {
                by_label[j] += c.by_label[j];
            }
        }
    }
};

struct HCountTable {

    const vector<vector<int> > &m_data;
    const vector<size_t> &m_labels;

    std::hash<Key> m_hasher;
    std::equal_to<Key> m_key_equal;

    // m_table must be declared after m_hasher and m_key_equal:
    std::unordered_map<Key, Count> m_table;

    HCountTable(
            const vector<vector<int> > &data,

            // The input labels must be consecutive starting from 0.
            //   E.g., 0, 1, 2, 3, ...
            const vector<size_t> &labels = vector<size_t>(0)
    )
            : m_data(data), m_labels(labels), m_hasher(data),
              m_key_equal(data), m_table(0, m_hasher, m_key_equal) {
        auto n = data.size();


        size_t nlabels = 0;

        if (labels.size() == n) {
            nlabels = 1 + *std::max_element(labels.begin(), labels.end());
        }

        if (nlabels == 0) {

            for (auto i = 0u; i < n; ++i) {

                Key key(i);

                auto p = m_table.find(key);

                if (p == m_table.end()) {

                    Count c(1);
                    c.members.push_back(i);

                    m_table.insert({key, c});
                } else {
                    p->second.tot++;
                    p->second.members.push_back(i);
                }
            }

        } else {

            for (auto i = 0u; i < n; ++i) {

                Key key(i);

                auto p = m_table.find(key);

                if (p == m_table.end()) {

                    Count c(1, nlabels);
                    c.by_label[labels[i]] = 1;
                    c.members.push_back(i);

                    m_table.insert({key, c});

                } else {
                    p->second.tot++;
                    p->second.by_label[labels[i]]++;
                    p->second.members.push_back(i);
                }
            }

        }
    }

    Count operator()(size_t i) noexcept {
        return m_table[Key(i)];
    }

    void print() {
        for (auto p = m_table.begin(); p != m_table.end(); ++p) {
            auto &v = m_data[(p->first).m_i];
            for (auto &x: v) cout << x << ",";
            cout << ": " << p->second.tot << ";\t";
            auto &l = p->second.by_label;
            for (auto &x: l) cout << x << ",";
            cout<< ";\t";
            auto &m = p->second.members;
            for (auto &x: m) cout << x << ",";
            cout << endl;
        }
    }
};

void test_HCountTable();

#endif /* HCountTable_h */
