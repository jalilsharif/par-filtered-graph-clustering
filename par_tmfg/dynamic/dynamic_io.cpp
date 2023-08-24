#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

#include "tmfg_dynamic.hpp"

using parlay::sequence;	




template <class T>
void dynamic_TMFG<T>::readValuesFromFile() {
    //timer t; t.start();
    std::ifstream myFile (fname, ios_base::in | ios_base::binary);
    for(int el=0; el<n; el++){
        myFile >> true_labels[el];
        max_label = std::max(max_label, true_labels[el]);
        for(int i=0; i<len; i++){
            myFile >> raw_data[el*len+i];

        }
    }

}


namespace IO_dynamic {
  using namespace std;

  extern ostream *time_output;

  template <class T, class Seq>
  SymM<T> parseSymMatrix(Seq W, std::size_t n) {
    SymM<T> matrix = SymM<T>(n);
    parlay::parallel_for(0, n, [&](size_t i){
        parlay::parallel_for(i+1, n,[&] (size_t j){
          matrix.update(i, j, (T)W[i*n + j]);
        });
    });
    return matrix;
  }

  // read a symmatric matrix from file


}  // namespace IO