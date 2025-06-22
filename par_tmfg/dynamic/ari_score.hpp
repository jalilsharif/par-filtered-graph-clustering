#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
 
using parlay::sequence;

template<class T> double ari_score(sequence<T> seq1, sequence<T> seq2, int num_clusters){
    int length = seq1.size();
    // to allow for 0-indexing or 1-indexing
    int slots = num_clusters + 1;

    sequence<int> pair_counter = sequence<int>(slots*slots, 0);
    sequence<int> seq1_counter = sequence<int>(slots, 0);
    sequence<int> seq2_counter = sequence<int>(slots, 0);
    for(int idx=0; idx<length; idx++){
        int seq1_label = seq1[idx];
        int seq2_label = seq2[idx];
        seq1_counter[seq1_label]++;
        seq2_counter[seq2_label]++;
        pair_counter[seq1_label*slots+seq2_label]++;
    }

    double sumpairs = 0;
    for(int i=0; i<slots*slots;i++){
        sumpairs += (pair_counter[i]*(pair_counter[i]-1))/2;
    }

    double sum1s = 0;
    for(int i=0; i<slots;i++){
        sum1s += (seq1_counter[i]*(seq1_counter[i]-1))/2;
    }

    double sum2s = 0;
    for(int i=0; i<slots;i++){
        sum2s += (seq2_counter[i]*(seq2_counter[i]-1))/2;
    }

    double prodsums = sum1s * sum2s;
    double sums = sum1s + sum2s;

    double npairs = length * (length - 1) / 2;

    double numerator = sumpairs - prodsums / npairs;
    double denominator = sums / 2 - prodsums / npairs;
    
    return numerator / denominator;
    

}