//
// Created by Jiandong Wang on 2022/11/16.
//

#ifndef GOCTEST_MYHASH_H
#define GOCTEST_MYHASH_H

#include "HCountTable.h"

struct MyHash : HCountTable {

    MyHash(const vector<vector<int> > &data,
            // The input labels must be consecutive starting from 0.
            //   E.g., 0, 1, 2, 3, ...
           const vector<size_t> &labels = vector<size_t>(0)) : HCountTable(data, labels) {}

    std::vector<vector<size_t>> gen_table() {
        size_t nlabels = 1 + *std::max_element(m_labels.begin(), m_labels.end());
        std::vector<vector<size_t>> conti_table(m_table.size(), std::vector<size_t>(nlabels, 0));
        auto iter_conti = conti_table.begin();

        //cluster is row, first loop
        for (auto p = m_table.begin(); p != m_table.end(); ++p) {
            *iter_conti = p->second.by_label;
            iter_conti++;
        }

/*        //print out for testing
        for(auto i:conti_table){
            for(auto j:i){
                cout<<j<<",";
            }
            cout<<endl;
        }*/
        return conti_table;
    }

};

#endif //GOCTEST_MYHASH_H
