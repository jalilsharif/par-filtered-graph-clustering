#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "parlay/parallel.h"
#include "parlay/primitives.h"

#include "symmatrix.h"


namespace IO {
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
  template <class T>
  SymM<T> readSymMatrixFromFile(char const *fname, std::size_t n) {
    //timer t; t.start();
    parlay::sequence<double> W = parlay::sequence<double>(n*n);
    std::ifstream myFile (fname, ios::in | ios::binary);
    myFile.read((char*)W.data(), sizeof(double) * n*n);
    if (W.size() == 0) {
      cout << "readPointsFromFile empty file" << endl;
      abort();
    }
    if (W.size() % n != 0) {
      cout << "readPointsFromFile wrong file type or wrong dimension" << endl;
      abort();
    }
    //cout<<t.next()<<'\n';
    auto z = parseSymMatrix<T>(W.cut(0,W.size()), n);
    //cout<<t.next()<<"aaa\n";
    return z;
  }

}  // namespace IO