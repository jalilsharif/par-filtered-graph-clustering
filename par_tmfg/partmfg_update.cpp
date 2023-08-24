#pragma once

#include <immintrin.h>

#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "parlay/delayed_sequence.h"
#include "parlay/hash_table.h"

#include "partmfg.h"

// Returns the vertex with max correlation to vertex i that hasn't been inserted yet into the TMFG
// Returns n+1 if there is no such vertex
// Has support for AVX512/AVX2
#ifdef __AVX512F__
template<class T, class PROF>
vtx ParTMFG<T, PROF>::getMaxCorr(vtx i){
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
        return n+1;
    }
    return ele;
}
#elif defined(__AVX2__)
template<class T, class PROF>
vtx ParTMFG<T, PROF>::getMaxCorr(vtx i){
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
        return n+1;
    }
    return ele;
}
#else
template<class T, class PROF>
vtx ParTMFG<T, PROF>::getMaxCorr(vtx i){
    vtx ele;
    size_t ind = 0;
    do {
        ind = corr_sorted_list_pointer[i];
        corr_sorted_list_pointer[i]++;
        ele = corr_list[i*n+ind];

    }
    while ( !vertex_flag[ele] && ind < n);
    if(ind == n && !vertex_flag[ele]){
        return n+1;
    }
    return ele;
}
#endif // AVX


template<class T, class PROF>
typename ParTMFG<T, PROF>::heapEle ParTMFG<T, PROF>::getApproxMaxGain(triT triangle, int depth){
    auto gainList = sequence<heapEle>::uninitialized(3 * depth);
    int counter = 0;
    
    vtx t1,t2,t3;
    tie(t1,t2,t3) = triangle;
    vtx vertices[3] = {t1, t2, t3};
    for(vtx v : vertices){
        size_t temp_ptr = corr_sorted_list_pointer[v];            
        for(int i = 0; i < depth; i++){
            vtx cur_elt = getMaxCorr(v);
            if(cur_elt == n+1){
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

// parallelizable, assumes gains are already calculated
template<class T, class PROF>
typename ParTMFG<T, PROF>::heapEle ParTMFG<T, PROF>::getFastMaxGain(triT triangle, int depth){
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

template<class T, class PROF> 
void ParTMFG<T, PROF>::updateVerticesHeap(vtx v, face tri){
    
    // Not parallelizable yet
    int depth = 1;

    pop_heap(vtx_heap.begin(), vtx_heap.end(), heap_compare());
    vtx_heap.pop_back();

    for(size_t ii = 0; ii < vtx_to_face_inds[v]; ii++){
        face i = face_store[ii+max_face_num*v];
        invalid_flag[i] = true;
        if(i == tri || i >= triangles_ind -2 ){
            triT new_tri = triangles[i];
            heapEle result = getApproxMaxGain(new_tri, depth);
            vtx new_v = result.second;
            max_clique_gains[i] = make_tuple(new_v, result.first, i);
            vtx_heap.push_back(max_clique_gains[i]);
            push_heap(vtx_heap.begin(), vtx_heap.end(), heap_compare());
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
            face_store[num_faces-1+max_face_num*new_v] = i;
            invalid_flag[i] = false;
        }
    }
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::updateVerticesCorrs(vtx v, face tri){
    int depth = 1;
    // get a set of all vertices to be updated
    auto update_set = parlay::hashtable(3 * vtx_to_face_inds[v], parlay::hash_numeric<T>());
    parlay::parallel_for(0, vtx_to_face_inds[v], [&](size_t ii) {
        face i = face_store[ii+max_face_num*v];
        triT new_tri = triangles[i];
        vtx v1, v2, v3;
        tie(v1,v2,v3) = new_tri;
        update_set.insert(v1);
        update_set.insert(v2);
        update_set.insert(v3);

    });
    // update the correlations in parallel
    sequence<T> vertices = update_set.entries();
    parlay::parallel_for(0, vertices.size(), [&](size_t vt) {
        getMaxCorr(vertices[vt]);
        corr_sorted_list_pointer[vertices[vt]]--;
    });
    // find the new gains in parallel
    parlay::parallel_for(0, vtx_to_face_inds[v], [&](size_t ii){
        face i = face_store[ii+max_face_num*v];
        triT new_tri = triangles[i];
        heapEle result = getFastMaxGain(new_tri, depth);
        vtx new_v = result.second;
        max_clique_gains[i] = make_tuple(new_v, result.first, i);
        max_gains_array[i] = (float)(result.first);
        size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
        face_store[num_faces-1+max_face_num*new_v] = i;
    });
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::updateVerticesOld(vtx v, face tri){
    parlay::parallel_for(0, vtx_to_face_inds[v], [&](size_t ii) {
            face i = face_store[ii+max_face_num*v];
            if(i == tri || i >= triangles_ind -2 ){heapifyFace(i);}

            heapEle result = getMinValidHeapEle(i);
            vtx new_v = result.second;
            max_clique_gains[i] = make_tuple(new_v, result.first, i);
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
            face_store[num_faces-1+max_face_num*new_v] = i;
    });
}