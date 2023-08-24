#pragma once

#include <tuple>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <algorithm>


#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/delayed_sequence.h"
#include "parlay/hash_table.h"


#ifdef HIGHWAY_MAKE
#include "hwy/base.h"
#include "hwy/contrib/sort/algo-inl.h"
#include "hwy/contrib/sort/result-inl.h"
#include "hwy/contrib/sort/traits128-inl.h"
#include <hwy/contrib/sort/vqsort.h>
#include <hwy/aligned_allocator.h>
#include <hwy/highway.h>
#include <hwy/nanobenchmark.h>
#include "hwy/tests/test_util-inl.h"
#endif // HIGHWAY_MAKE


#include "gettime.h"
#include "atomics.h"
#include "utilities.h"
#include "IO.h"

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

    ParTMFG(SymM<T> *W_, size_t n_, PROF *_profiler, bool _use_corrs, bool _use_gains_heap, bool _use_highway):
    W(W_), n(n_), pf(_profiler), use_corrs(_use_corrs), use_max_gains_heap(_use_gains_heap), use_highway(_use_highway){
    }

    struct heap_compare{
    bool operator()(const gainT& i,const gainT& j) const{
        return get<1>(i) < get<1>(j);
        }
    };




    // initialize the four cliques
    void init();


    // start of insertOne
    gainT getMaxVtxArray();
    gainT getMaxVtxHeap();

    void updateInternal(vtx v, face tri, DBHTTMFG<T, PROF> *clusterer);

    // functions at end of insertOne
    void updateVerticesHeap(vtx v, face tri);
    void updateVerticesCorrs(vtx v, face tri);
    void updateVerticesOld(vtx v, face tri);

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
    void insertMultiple(sequence<size_t> &insert_list, DBHTTMFG<T, PROF> *clusterer=nullptr);

    T computeGain(vtx i, triT triangle){
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

    void insertToP(vtx i, vtx j, size_t ind){
        P[ind] = make_tuple(i, j, getW(i, j));
    }

    bool hasUninsertedV(){return vertex_num > 0;}  

    T getW(vtx i, vtx j){
        return W->get(i,j);
    }

    size_t getTrianglesNum(){return triangles_ind;}

    T computeCost(){
      // compute the total gain here
        T total_cost = 0;
        for(size_t i=0; i<P_ind; ++ i){
            total_cost += get<2>(P[i]);
        }

        (*IO::time_output) << std::setprecision(20) << 2*total_cost << endl;
        (*IO::time_output) << P_ind << endl; 
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
    
    sequence<heapEle> heap_buffer;
    sequence<heapEle> corr_heap_buffer;

    sequence<vtx> corr_list;

    sequence<size_t> sorted_list_pointer; // stores the next vertex to look at
    sequence<size_t> corr_sorted_list_pointer;

    vector<gainT> vtx_heap;
    sequence<bool> invalid_flag;
    bool use_corrs;
    bool use_max_gains_heap;
    bool use_highway;


    sequence<float> max_gains_array;



    //allocate space for heap
    void initGainArrayHeap();
    void updateGainArrayHeap(sequence<size_t> &insert_list);

    heapEle negateGain(heapEle ele){
        return heapEle(-1.0*ele.first, ele.second);
    }

    //prereq: vertex_flag is updated
    inline heapEle getMinValidHeapEle(face i){
        pf->incFindMinNum();
        heapEle ele;
        do {
            size_t ind = sorted_list_pointer[i];
            sorted_list_pointer[i]++;
            ele = heap_buffer[i*n+ind];
        }
        while ( !vertex_flag[ele.second] );
        return negateGain(ele);
    }


    vtx getMaxCorr(vtx i);


    // finds the closest uninserted vertex to each of the face's three vertices
    // and calculates the uninserted vertex with the max gain
    heapEle getApproxMaxGain(triT triangle, int depth);

    // parallelizable
    heapEle getFastMaxGain(triT triangle, int depth);

    
    // init a heap buffer for triangles[i]
    void heapifyFace(face i);

    // sorts the vertices by correlation with a given vertex
    void heapifyVtx(vtx i);



    //// used for heap

    //tested
    cliqueT maxClique();

    void removeOneV(vtx v);
    void filterVbyFlag();

};

