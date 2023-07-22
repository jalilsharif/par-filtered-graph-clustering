#pragma once

#include <tuple>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <algorithm>
#include <immintrin.h>


#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/delayed_sequence.h"
#include "parlay/hash_table.h"

#include <boost/sort/spreadsort/float_sort.hpp>

#define VQSORT_ONLY_STATIC 1
#define VQSORT_PRINT 0

#include "hwy/base.h"

#include "hwy/contrib/sort/algo-inl.h"
#include "hwy/contrib/sort/result-inl.h"
#include "hwy/contrib/sort/traits128-inl.h"
#include <hwy/contrib/sort/vqsort.h>
#include <hwy/aligned_allocator.h>
#include <hwy/highway.h>
#include <hwy/nanobenchmark.h>
#include "hwy/tests/test_util-inl.h"


#include "gettime.h"
#include "atomics.h"
#include "utilities.h"

#include "heap.h"
#include "dbht.h"

using namespace std;
using parlay::sequence;	
using parlay::slice;	
using parlay::hashtable;
using parlay::hash_numeric;

template<class T, class PROF> 
struct ParTMFG{
    using gainT=tuple<vtx, T, face>; 
    using heapEle=pair<T, vtx>;
    using heapT=binary_min_heap<heapEle>;
    // struct entry {
    // using key_t = heapEle;
    // static bool comp(key_t a, key_t b) { 
    //     return a < b;}
    // };
    // using heapT=pam_set<entry>;
    SymM<T>* W;
    sequence<tuple<vtx,vtx,T>> P;
    sequence<gainT> max_clique_gains;

    sequence<cliqueT> cliques;
    sequence<triT> triangles;
    sequence<vtx> peo;
    size_t n;

    PROF *pf;

    ParTMFG(SymM<T> *W_, size_t n_, PROF *_profiler, bool _use_corrs, bool _use_gains_heap, bool _use_highway, bool _use_heap=false):
    W(W_), n(n_), pf(_profiler), use_corrs(_use_corrs), use_max_gains_heap(_use_gains_heap), use_highway(_use_highway), use_heap(_use_heap){
    }

    struct heap_compare{
    bool operator()(const gainT& i,const gainT& j) const{
        return get<1>(i) < get<1>(j);
        }
    };

    struct rightshift {
         long long operator()(heapEle e, unsigned offset) { return boost::sort::spreadsort::float_mem_cast<T, long long>(e.first * 8192) >> offset; }
    };



    // initialize the four cliques
    void init();

    //require max_clique_gain updated
    void insertOne(DBHTTMFG<T, PROF> *clusterer=nullptr);
    void initGainArray();


    // max_clique_gain does not neccessarily correspond to the order in triangles anymore
    parlay::sequence<size_t> getAllBestVertices(size_t THRESHOLD);

    // get the best vertices of the largest P gains in max_clique_gain vector
    // require: THRESHOLD must be not larger than triangles_ind
    // max_clique_gain does not neccessarily correspond to the order in triangles anymore
    parlay::sequence<size_t> getBestVertices(size_t THRESHOLD);

    // update the max_clique_gain array after calling 
    // update the faces in [0, vtx_to_face_inds[v]] of face_store[v*n] for all v in insert_list
    // max_clique_gain does not neccessarily correspond to the order in triangles when start, but does after
    void updateGainArray(sequence<size_t> &insert_list);

    // insert all vertices in max_clique_gain[insert_list].
    // i.e. the indices in insert_list are indices of max_clique_gain
    void inertMultiple(sequence<size_t> &insert_list, DBHTTMFG<T, PROF> *clusterer=nullptr);

    inline T computeGain(vtx i, triT triangle){
        return getW(i,get<0>(triangle)) + getW(i,get<1>(triangle)) + getW(i,get<2>(triangle));
    }

    //get the best vertex for triangle in vertex_list
    pair<vtx, T> getBestGain(parlay::slice<vtx *, vtx *> &vertex_list, triT triangle){
// #ifdef PROFILE
//         timer t;t.start();
// #endif  
        vtx *best_vertex = parlay::max_element(vertex_list, 
            [&](const auto &i, const auto &j){ 
                T ga = computeGain(i,triangle);
                T gb = computeGain(j,triangle); 
                if(ga == gb){
                    return i < j;
                    //want the smallest id vertex if gain equal in matlab, i > j
                }
                return  ga < gb;
            });
        T gain = computeGain(*best_vertex,triangle);
        return make_pair(*best_vertex, gain);
    }

    inline void insertToP(vtx i, vtx j, size_t ind){
        P[ind] = make_tuple(i, j, getW(i, j));
    }

    inline bool hasUninsertedV(){return vertex_num > 0;}    
    inline T getW(vtx i, vtx j){
// #ifdef DEBUG
//         if(i==j){
//             cout << "i==j in getW" << endl;
//         }
// #endif
        //{return W[i*n + j];}
        return W->get(i,j);
    }
    inline size_t getTrianglesNum(){return triangles_ind;}

    T computeCost(){
        ofstream outfile;
        outfile.open("timedata.txt", std::ios_base::app);
      // compute the total gain here
        T total_cost = 0;
        for(size_t i=0; i<P_ind; ++ i){
            // print(P[i]);
            total_cost += get<2>(P[i]);
        }

        outfile << std::setprecision(20) << 2*total_cost << endl;
        outfile << P_ind << endl; 
        return 2*total_cost;   
    }

    void outputPeo(string filename, size_t _n=0){
        if(_n==0)_n = n;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << peo[i]+1 << endl;
        }
        file_obj.close();
    }    

    void outputP(string filename, size_t _n=0){
        if(_n==0)_n = P_ind;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << get<0>(P[i])+1 << " " << get<1>(P[i])+1 << " " << get<2>(P[i]) << endl;
        }
        file_obj << n << " " << n << " " << 0 << endl;
        file_obj.close();
    }   

    void outputCliques(string filename, size_t _n=0){
        if(_n==0)_n = n-3;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << get<0>(cliques[i]) << " " << get<1>(cliques[i]) << " " << get<2>(cliques[i]) << " " << get<3>(cliques[i]) << endl;
        }
        file_obj.close();
    }    
#ifdef DEBUG
        bool validV(vtx t){
            return t < n && t >= 0;
        }

        void checkTriangles(){
            for (size_t i=0; i < triangles_ind; ++i){
                vtx t1,t2,t3;
                tie(t1,t2,t3) = triangles[i];
                if(t1 == t2 || t1==t3 || t2==t3){
                    cout << t1 << " " << t2 << " " << t3 << endl;
                }
            }


        }

        void checkGains(size_t k){
            for (size_t i=0; i < triangles_ind-k; ++i){
                auto entry = max_clique_gains[i];
                face tri = get<2>(entry);
                if(tri != i){
                    cout << tri << " " << i << endl;
                }
            }


        }
#endif

    private:
    parlay::sequence<bool> vertex_flag;
    sequence<vtx> vertex_list;
    sequence<size_t> face_store;
    // sequence<T> W2;

    sequence<size_t> vtx_to_face_inds; //how many faces have vtx as the best vtx
    size_t vertex_num;
    size_t vertex_start;
    size_t triangles_ind;
    size_t peo_ind;
    size_t P_ind = 0;
    size_t max_face_num;
    //// used for heap
    bool use_heap;
    sequence<heapT> heaps;
    sequence<heapEle> heap_buffer;
    sequence<heapEle> corr_heap_buffer;

    sequence<vtx> corr_list;

    sequence<size_t> heap_LR;
    bool use_sorted_list = true; // change to input
    sequence<size_t> sorted_list_pointer; // stores the next vertex to look at
    sequence<size_t> corr_sorted_list_pointer;

    vector<gainT> vtx_heap;
    sequence<bool> invalid_flag;
    bool use_corrs;
    bool use_max_gains_heap;
    bool use_highway;

    sequence<size_t> corr_second_pointer;

    sequence<float> max_gains_array;
    sequence<vtx> max_vtx_array;

    

    size_t block_freq = 500;
    size_t add_buffer = 5;




    //allocate space for heap
    void initHeap();
    void initGainArrayHeap();
    void updateGainArrayHeap(sequence<size_t> &insert_list);

    inline heapEle negateGain(heapEle ele){
        return heapEle(-1.0*ele.first, ele.second);
    }

    //prereq: vertex_falg is updated
    inline heapEle getMinValidHeapEle(face i){
// #ifdef PROFILE
//         timer t;t.start();
// #endif  
        pf->incFindMinNum();
        heapEle ele;
        do {
// #ifdef DEBUG
//     if(heaps[i].size == 0){
//         cout << "0 size heap!" <<endl;
//     }
// #endif
            if(use_sorted_list){ 
                size_t ind = sorted_list_pointer[i];
                sorted_list_pointer[i]++;
                ele = heap_buffer[i*n+ind];
            }else{
            size_t ind = heaps[i].argmin();
            ele = heap_buffer[i*n+ind];
            heaps[i].heap_pop();
            // ele = *heaps[i].select(0);
            // heaps[i]=heapT::remove(move(heaps[i]), ele);
            }
        }
        while ( !vertex_flag[ele.second] );
// #ifdef PROFILE
//         // pf->incFindMinTime(t.next());
//         cout << "find min: "<<t.next()<<endl;
// #endif  
        return negateGain(ele);
    }


    // Returns the vertex with max correlation to vertex i that hasn't been inserted yet into the TMFG
    // Returns -1 if there is no such vertex
    // Has support for AVX512/AVX2
    #if __AVX512F__
    inline vtx getMaxCorr(vtx i){
        vtx ele;
        size_t ind = 0;
        ind = corr_sorted_list_pointer[i];
        
        corr_sorted_list_pointer[i]++;
        ele = corr_list[i*n+ind];
        if(vertex_flag[ele]){
            return ele;
        }
        ind = corr_sorted_list_pointer[i];
        while(ind < n - 17){
            __m512i els = _mm512_loadu_si512((__m512i_u*)(corr_list.data() + i * n + ind));
            __m512i flags = _mm512_i32gather_epi32(els, (int*)(vertex_flag.data()), sizeof(bool));
            __mmask16 masked_flags = _mm512_cmpneq_epi32_mask((flags & _mm512_set1_epi32(1)), _mm512_set1_epi32(0));

            if((int)masked_flags != 0){
                int first_idx = __builtin_ctz((int)masked_flags);
                corr_sorted_list_pointer[i] += (first_idx + 1);
                return corr_list[i * n + ind + first_idx];
            }
            corr_sorted_list_pointer[i] += 16;
            ind = corr_sorted_list_pointer[i];
            
        }

        while(!vertex_flag[ele] && ind < n){
            ind = corr_sorted_list_pointer[i];
            corr_sorted_list_pointer[i]++;
            ele = corr_list[i*n+ind];

        }
        if(ind == n && !vertex_flag[ele]){
            return -1;
        }
        return ele;
    }
    #elif __AVX2__
    inline vtx getMaxCorr(vtx i){
        vtx ele;
        size_t ind = 0;
        ind = corr_sorted_list_pointer[i];
        
        corr_sorted_list_pointer[i]++;
        ele = corr_list[i*n+ind];
        if(vertex_flag[ele]){
            return ele;
        }
        ind = corr_sorted_list_pointer[i];
        while(ind < n - 9){
            __m256i els = _mm256_loadu_si256((__m256i_u*)(corr_list.data() + i * n + ind));
            __m256i flags = _mm256_i32gather_epi32((int*)(vertex_flag.data()), els, sizeof(bool));
            
            int masked_flags = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32((flags & _mm256_set1_epi32(1)), _mm256_set1_epi32(1))));

            
            if(masked_flags != 0){
                int first_idx = __builtin_ctz((int)masked_flags);
                corr_sorted_list_pointer[i] += (first_idx + 1);
                return corr_list[i * n + ind + first_idx];
            }
            corr_sorted_list_pointer[i] += 8;
            ind = corr_sorted_list_pointer[i];
            
        }

        while(!vertex_flag[ele] && ind < n){
            ind = corr_sorted_list_pointer[i];
            corr_sorted_list_pointer[i]++;
            ele = corr_list[i*n+ind];

        }
        if(ind == n && !vertex_flag[ele]){
            return -1;
        }
        return ele;
    }
    #else
    inline vtx getMaxCorr(vtx i){
        vtx ele;
        size_t ind = 0;
        do {
            ind = corr_sorted_list_pointer[i];
            corr_sorted_list_pointer[i]++;
            ele = corr_list[i*n+ind];

        }
        while ( !vertex_flag[ele] && ind < n);
        if(ind == n && !vertex_flag[ele]){
            return -1;
        }
        return ele;
    }
    #endif // AVX


    // finds the closest uninserted vertex to each of the face's three vertices
    // and calculates the uninserted vertex with the max gain
    inline heapEle getApproxMaxGain(triT triangle, size_t depth){
        auto gainList = sequence<heapEle>::uninitialized(3 * depth);
        int counter = 0;
        
        vtx t1,t2,t3;
        tie(t1,t2,t3) = triangle;
        vtx vertices[3] = {t1, t2, t3};
        for(vtx v : vertices){
            size_t temp_ptr = corr_sorted_list_pointer[v];            
            for(int i = 0; i < depth; i++){
                vtx cur_elt = getMaxCorr(v);
                if(cur_elt == -1){
                    break; 
                }
                if(i == 0){
                    temp_ptr = corr_sorted_list_pointer[v] - 1;
                }
                gainList[counter] = heapEle(computeGain(cur_elt, triangle), cur_elt);
                counter++;
                
            }
            corr_sorted_list_pointer[v] = temp_ptr;
        }
        return *parlay::max_element(make_slice(gainList).cut(0, counter), std::less<heapEle>{});
    }

    // parallelizable
    inline heapEle getFastMaxGain(triT triangle, size_t depth){
        auto gainList = sequence<heapEle>::uninitialized(3 * depth);
        int counter = 0;
        
        vtx t1,t2,t3;
        tie(t1,t2,t3) = triangle;
        vtx vertices[3] = {t1, t2, t3};
        for(vtx v : vertices){
            for(int i = 0; i < depth; i++){
                if(corr_sorted_list_pointer[v] >= n){
                    break;
                }
                vtx cur_elt = corr_list[v * n + corr_sorted_list_pointer[v]];
                gainList[counter] = heapEle(computeGain(cur_elt, triangle), cur_elt);
                counter++;
                
            }
        }
        return *parlay::max_element(make_slice(gainList).cut(0, counter), std::less<heapEle>{});
    }


    // compresses the correlation list to remove some uninserted vertices
    inline void reorganize_vtx(vtx v){
        int valids = 0;
        corr_sorted_list_pointer[v] = 0;
        for(int i = 0; i < block_freq + add_buffer; i++){
            vtx ele = corr_list[v*n+i];
            if(vertex_flag[ele]){
                corr_list[v*n+valids] = ele;
                valids++;
            }
        }
        for(int i = corr_second_pointer[v]; i<n; i++){
            vtx ele = corr_list[v*n+i];
            corr_second_pointer[v]++;
            if(vertex_flag[ele]){
                corr_list[v*n+valids] = ele;
                valids++;
                if(valids >= block_freq + add_buffer){
                    break;
                }
            }
            
        }
    }



    
    // init a heap buffer for triangles[i]
    inline void heapifyFace(face i){
// #ifdef PROFILE
//         timer t1;t1.start();
// #endif  
        pf->incHeapifyNum();
        pf->incHeapifySize(vertex_num);
        triT t = triangles[i];
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
        /*parlay::parallel_for(0, in.size(), [&](size_t v_ind) {
            vtx v = in[v_ind];
            T gain = computeGain(v, t);
            heap_buffer[i*n+v_ind] = heapEle(-1.0*gain, v);
        });*/
        if(use_sorted_list){
            if(use_highway){
                vector<uint64_t> z = vector<uint64_t>(vertex_num);

                for(size_t v = 0; v < vertex_num; v++){
                    T gain = -1.0 * computeGain(in[v], t);
                    //corr_heap_buffer[i*n+v] = heapEle(-1.0*gain, v);
                    float fl_gain = (float)(gain+3);
                    uint32_t f_gain;
                    memcpy(&f_gain, &(fl_gain), sizeof(float));
                    //f_gain+=v;
                    z[v] = (((uint64_t)f_gain) << 32) + in[v];
                }   
                // parlay::sort_inplace(make_slice(heap_buffer).cut(i*n, i*n+vertex_num));
                // parlay::internal::seq_sort_inplace(make_slice(heap_buffer).cut(i*n, i*n+vertex_num), std::less<heapEle>{}, false);
                //parlay::internal::quicksort(make_slice(heap_buffer).cut(i*n, i*n+vertex_num), std::less<heapEle>{});
                hwy::VQSort(&z[0], vertex_num, hwy::SortAscending());
                for(vtx v = 0; v < vertex_num; v++){
                    //float d_gain;
                    //uint32_t u_gain = (uint32_t)(z[v] >> 32);
                    //memcpy(&d_sgain, &(u_gain), 4);
                    uint32_t f_gain = (uint32_t)(z[v]>>32);
                    float fl_gain;
                    memcpy(&fl_gain, &(f_gain), sizeof(float));
                    heap_buffer[i*n+v]=heapEle((T)(fl_gain-3), (uint32_t)(z[v]));
                }
            }
            else{
                parlay::parallel_for(0, in.size(), [&](size_t v_ind) {
                    vtx v = in[v_ind];
                    T gain = computeGain(v, t);
                    heap_buffer[i*n+v_ind] = heapEle(-1.0*gain, v);
                });
                parlay::internal::quicksort(make_slice(heap_buffer).cut(i*n, i*n+vertex_num), std::less<heapEle>{});
            }

            sorted_list_pointer[i] = 0;

       }else{
        heaps[i] = binary_min_heap<heapEle>(heap_buffer.data()+ (i*n), vertex_num, heap_LR.data()+ (i*n));
        heaps[i].heapify();
        // heaps[i] = heapT();
        // heaps[i] = heapT::multi_insert(move(heaps[i]), make_slice(heap_buffer).cut(i*n, i*n+vertex_num));
       }

// #ifdef PROFILE
//         // pf->incHeapifyTime(t1.next());
//         cout << vertex_num << endl;
//         cout << "heapify: "<<t1.next() << endl;;
// #endif  


    }

    // sorts the vertices by correlation with a given vertex
    inline void heapifyVtx(vtx i){
        vector<uint64_t> z = vector<uint64_t>(n);
        for(vtx v = 0; v < n; v++){
            T corr = getW(i, v);
            // since correlations are between -1 and 1, add 1 to ensure corr is non-negative 
            // so conversion to unsigned integer preserves order
            float fl_corr = (float)(-1.0*corr+1);
            uint32_t f_corr;
            memcpy(&f_corr, &(fl_corr), sizeof(float));
            z[v] = (((uint64_t)f_corr) << 32) + v;
        }
        if(use_highway){
            hwy::VQSort(&z[0], n, hwy::SortAscending());
        }
        else{
            std::sort(&z[0], &z[0] + n);
        }

        corr_sorted_list_pointer[i] = 0;
        for(vtx v = 0; v < n; v++){
            // chop off the correlation weight to leave only the vertex number
            corr_list[i*n+v]=(uint32_t)(z[v]);
        }
    }



    //// used for heap

    //tested
    cliqueT maxClique();

    void removeOneV(vtx v){
        //update vertex list
        //filter and then copy back
        size_t new_vertex_start = 0;
        if(vertex_start == 0){
            new_vertex_start = n;
        }
        auto out2 = make_slice(vertex_list).cut(new_vertex_start, n);
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);

        vertex_num = parlay::filter_into(in, out2, [&](vtx i){ 
            return i != v; 
            });
        vertex_start = new_vertex_start;
    }

    void filterVbyFlag(){
        //update vertex list
        //filter and then copy back
        size_t new_vertex_start = 0;
        if(vertex_start == 0){
            new_vertex_start = n;
        }
        auto out2 = make_slice(vertex_list).cut(new_vertex_start, n);
        auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);

        vertex_num = parlay::filter_into(in, out2, [&](vtx i){ 
            return vertex_flag[i]; 
            });
        vertex_start = new_vertex_start;
    }

    /*void compareParent(face i, face tri){
        vtx best = getMinValidHeapEle(i).second;
        int k = 0;
        size_t ind = 0;

        heapEle ele;
        do {

                ele = heap_buffer[tri*n+ind];
                ind++;
                if(vertex_flag[ele.second]){
                    k++;
                    
                }
                if(ele.second==best){
                        break;
                    }

        }
        while ( ind<n );
        cout<<k<<' ';
    }*/

    /*void bestIdx(vtx start, vtx best){
        int k = 0;
        size_t ind = 0;

        heapEle ele;
        do {

                ele = corr_heap_buffer[start*n+ind];
                ind++;
                if(vertex_flag[ele.second]){
                    k++;
                    
                }
                if(ele.second==best){
                        break;
                    }

        }
        while ( ind<n );
        cout<<k<<' ';
    }*/

    //https://en.algorithmica.org/hpc/algorithms/argmin/
    int argmin2(int *a, int n) {

        if(n<32){
            int k = 0;

            for (int i = 0; i < n; i++)
                if (a[i] > a[k])
                    k = i;

            return k;
        }

        int max = INT_MIN, idx = 0;
        
        
        __m256i p = _mm256_set1_epi32(max);

        for (int i = 0; i < n - (n % 32); i += 32) {
            __m256i y1 = _mm256_load_si256((__m256i*) &a[i]);
            __m256i y2 = _mm256_load_si256((__m256i*) &a[i + 8]);
            __m256i y3 = _mm256_load_si256((__m256i*) &a[i + 16]);
            __m256i y4 = _mm256_load_si256((__m256i*) &a[i + 24]);
            y1 = _mm256_max_epi32(y1, y2);
            y3 = _mm256_max_epi32(y3, y4);
            y1 = _mm256_max_epi32(y1, y3);
            __m256i mask = _mm256_cmpgt_epi32(y1, p);
            if (!_mm256_testz_si256(mask, mask)) { [[unlikely]]
                idx = i;
                for (int j = i; j < i + 32; j++)
                    max = (a[j] > max ? a[j] : max);
                p = _mm256_set1_epi32(max);
            }
        }

        int k = idx + 31;

        


        for (int i = idx; i < idx + 32; i++)
            if (a[i] == max){
                k = i;
                break;
            }

        for (int i = n - (n % 32); i < n; i++)
            if (a[i] > a[k])
                k = i;

        return k;
    }
};

