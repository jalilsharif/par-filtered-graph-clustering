#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

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
            parlay::parallel_for(i+1, n, [&](size_t j){
                matrix.update(i, j, (T)W[i*n + j]);
            });
        });
        return matrix;
    }

    template <class T>
    SymM<T> readSymMatrixFromFile(const char* fname, std::size_t n) {
        try {
            // Validate input arguments
            if (!fname || n == 0) {
                throw invalid_argument("Invalid file name or matrix size.");
            }

            // Open the file in binary mode
            ifstream myFile(fname, ios::in | ios::binary);
            if (!myFile.is_open()) {
                throw runtime_error("Error: Cannot open file " + string(fname));
            }

            // Validate file size
            myFile.seekg(0, ios::end);
            size_t file_size = myFile.tellg();
            myFile.seekg(0, ios::beg);

            size_t expected_size = sizeof(double) * n * n;
            if (file_size != expected_size) {
                throw runtime_error("File size mismatch: File size = " + to_string(file_size) +
                                    ", Expected size = " + to_string(expected_size));
            }

            // Load the binary data into memory
            parlay::sequence<double> W(n * n);
            myFile.read(reinterpret_cast<char*>(W.data()), expected_size);

            if (!myFile.good()) {
                throw runtime_error("Error reading file data from " + string(fname));
            }

            // Close the file
            myFile.close();

            // Debugging output for verification
            cout << "Successfully loaded matrix of size: " << n << " x " << n << endl;
            cout << "First few elements of the matrix: ";
            for (size_t i = 0; i < std::min((size_t)10, W.size()); ++i) {
                cout << W[i] << " ";
            }
            cout << endl;

            return parseSymMatrix<T>(W.cut(0, W.size()), n);

        } catch (const exception& e) {
            cerr << "Exception in readSymMatrixFromFile: " << e.what() << endl;
            exit(1); // Fail fast on critical error
        } catch (...) {
            cerr << "Unknown error occurred in readSymMatrixFromFile." << endl;
            exit(1); // Fail fast on unexpected error
        }
    }
} // namespace IO
