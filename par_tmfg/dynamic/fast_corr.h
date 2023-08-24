#pragma once
#include <iostream>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

#include "../symmatrix.h"


using parlay::sequence;


//div by 0 in corr calculations??????????

template<typename T>
struct fast_timeseries{
    int window_size;
    int min_idx;
    //int timeseries_size;
    //T* data;
    sequence<T> data;
    T sliding_sum;
    T sliding_square_sum;

    fast_timeseries(T* _data, int _window_size, int start, int len):
    window_size(_window_size), min_idx(start){
        data = sequence<T>::uninitialized(len);
        for(int i=0; i<len;i++){
            data[i]= _data[i];
        }
        sliding_sum = 0; 
        sliding_square_sum = 0;

        for(int i = min_idx; i < min_idx + window_size; i++){
            sliding_sum += data[i];
            //std::cout<<data[i]<<'\n';
            sliding_square_sum += data[i] * data[i];
        }
        //std::cout<<sliding_square_sum<<"!\n";
    }
    
    void advance();
};


template<typename T>
struct fast_timeseries_list{
    sequence<fast_timeseries<T>> series;
    int num_series;
    int len;
    sequence<T> sum_prods;
    int window_size;
    int min_idx;
    SymM<T> corr_matrix;
    

    fast_timeseries_list(){}

    fast_timeseries_list(T* raw_series, int _num_series, int _len, int _window_size, int start):
    num_series(_num_series), len(_len), window_size(_window_size), min_idx(start){
        series = sequence<fast_timeseries<T>>::uninitialized(num_series);

        for(int i = 0; i < num_series; i++){
            series[i] = fast_timeseries<T>(raw_series+i*len, window_size, min_idx, len);
            //sum_prods[i] = (T*) malloc(num_series * sizeof(T));
            //corr_matrix[i] = (double*) malloc(num_series * sizeof(double));
        }

        std::cout<<'\n';

    }

    void init(){
        sum_prods = sequence<T>::uninitialized(num_series * num_series);
        corr_matrix = SymM<double>(num_series);
        corr_matrix.setDiag(1);
        
        parlay::parallel_for(0, num_series, [&](int i){
            for(int j = 0; j < num_series; j++){
                corr_matrix.update(i,j,0);
                sum_prods[i*num_series+j] = 0;
                
                for(int k = min_idx; k < min_idx + window_size; k++){
                    sum_prods[i*num_series+j] += series[i].data[k] * series[j].data[k];
                    //std::cout<< series[i].data[k]<<' ';
                }
            }
        });

    }

    void advance();



    double compute_corr(int i, int j);

    void update_corr_matrix(){
        parlay::parallel_for(0, num_series, [&](int i){
            for(int j = i+1; j < num_series; j++){
                corr_matrix.update(i, j, compute_corr(i, j));

                //std::cout << corr_matrix.get(i,j)<<' ';
            }
            //std::cout<<'\n';
        });
        //std::cout<<'\n';
    }
        
};