#pragma once

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

#include "fast_corr.cpp"

#include <fstream>

#include "../symmatrix.h"

using parlay::sequence;


template<class T>
struct dynamic_TMFG{

    


    unsigned int n;
    unsigned int len;
    string fname1;
    string fname2;
    string fname;

    int time = 0;
    int window;

    int max_label = 0;

    int recluster_ct = 0;

    dynamic_TMFG(unsigned int n_, unsigned int len_, int window_, bool two_files, string fname1_, string fname2_, string fname_): 
        n(n_), len(len_), window(window_), fname1(fname1_), fname2(fname2_), fname(fname_){
        raw_data = sequence<T>::uninitialized(n*len);
        true_labels = sequence<int>::uninitialized(n);
        if(two_files){
            std::ifstream ifile1(fname1.c_str(), ios_base::in | ios_base::binary);
            std::ifstream ifile2(fname2.c_str(), ios_base::in | ios_base::binary);
    
            std::ofstream ofile(fname.c_str(), std::ios::out | ios_base::binary);
            ofile << ifile1.rdbuf() << ifile2.rdbuf();
            ofile.close();

        }
        readValuesFromFile();
        corrs = fast_timeseries_list<double>((T*)(raw_data.data()), n, len, window, 0); //(double*)(raw_data.data())
        corrs.init();
        prev_labels = sequence<int>(n);
    }



    void readValuesFromFile();

    void tick(bool start, bool stop, int i, int clusters, bool verbose, string method, int THRESHOLD, bool use_corrs, bool use_gains_heap, bool use_highway, bool exact_apsp, string dsname);

    sequence<int> runDBHT(SymM<double> *W, SymM<double> *D, size_t n, size_t THRESHOLD, string method, bool use_corrs, bool use_gains_heap, bool use_highway, bool exact_apsp, int num_clusters, int time, bool relabel, sequence<int> prev_labels, bool verbose, string dsname = "");

    T getVal(int el, int i){
        return raw_data[el*len + i];
    }

    private:
        sequence<int> true_labels;
        sequence<T> raw_data;
        fast_timeseries_list<double> corrs;
        sequence<int> prev_labels;


};