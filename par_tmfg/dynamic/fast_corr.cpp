#pragma once

#include "fast_corr.h"
#include <math.h>

template<typename T>
void fast_timeseries<T>::advance(){
    sliding_sum += data[min_idx + window_size];
    sliding_sum -= data[min_idx];
    sliding_square_sum += data[min_idx + window_size] * data[min_idx + window_size];
    sliding_square_sum -= data[min_idx] * data[min_idx];
    min_idx++;
}

template<typename T>
void fast_timeseries_list<T>::advance(){
    for(int i = 0; i < num_series; i++){
        series[i].advance();
    }
    parlay::parallel_for(0, num_series, [&](int i){
        for(int j = 0; j < num_series; j++){
            sum_prods[i*num_series+j] += series[i].data[min_idx + window_size] * series[j].data[min_idx + window_size];
            sum_prods[i*num_series+j] -= series[i].data[min_idx] * series[j].data[min_idx];
        }
    });
    //cout<<min_idx<<"aaaaaaaa\n";
    min_idx++;
}

template<typename T>
double fast_timeseries_list<T>::compute_corr(int i, int j){
    /*if(i==0&&j==2){
        std::cout<<sum_prods[i][j]<<' '<<series[i].sliding_sum<<' ' <<series[j].sliding_sum<<'\n';
    }*/
    double numerator = (double) (window_size * sum_prods[i*num_series+j] - series[i].sliding_sum * series[j].sliding_sum);
    double sqrt_i = sqrt((double) (window_size * series[i].sliding_square_sum - series[i].sliding_sum * series[i].sliding_sum));
    double sqrt_j = sqrt((double) (window_size * series[j].sliding_square_sum - series[j].sliding_sum * series[j].sliding_sum));
    /*std::cout<<(series[i].sliding_sum * series[j].sliding_sum)<<"?\n";
    std::cout<<(window_size * sum_prods[i*num_series+j])<<"?\n";
    std::cout<<sqrt_i<<"?\n";
    std::cout<<sqrt_j<<"?\n";*/
    if((double) (window_size * series[j].sliding_square_sum - series[j].sliding_sum * series[j].sliding_sum)<=0){
        cout<<"failure :("<<min_idx<<' '<<j;
    }
    return numerator / (sqrt_i * sqrt_j);
}

/*int main(){
    int window_size = 3;
    int num_series = 3;
    int start = 0;

    int arr_data[15] = {1, 2, 3, 4, 5, 
                        7, 6, 3, 4, 5, 
                        1, 2, 3, -8, 5};

    //std::cout<<z[0][0];
    //std::cout<<data[0][0]+data[0][1]+data[0][2]<<'\n';
    struct fast_timeseries_list<int> test_list = fast_timeseries_list<int>((int*)arr_data, num_series, 5, window_size, start);

    test_list.update_corr_matrix();
    test_list.advance();
    test_list.update_corr_matrix();
    test_list.advance();
    test_list.update_corr_matrix();

}*/