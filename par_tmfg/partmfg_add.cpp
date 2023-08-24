#pragma once

#include "parlay/primitives.h"
#include "parlay/sequence.h"

#include "partmfg.h"

// get the best vertex among all candidates using max_clique_gains
template<class T, class PROF> 
typename ParTMFG<T, PROF>::gainT ParTMFG<T, PROF>::getMaxVtxArray(){
    auto *entry_pointer = parlay::max_element(make_slice(max_clique_gains).cut(0,triangles_ind), 
                    [&](const auto &i, const auto &j){ 
                        if(get<1>(i) == get<1>(j)) return i > j;
                        // might not be exactly the same as MATLAB, 
                        // because triangle id can be different
                        return get<1>(i) < get<1>(j); 
                    });
    return *entry_pointer;
}

// get the best vertex among all candidates using heap
template<class T, class PROF> 
typename ParTMFG<T, PROF>::gainT ParTMFG<T, PROF>::getMaxVtxHeap(){
    int depth = 1;

    gainT entry;
    while(true){
        entry = vtx_heap.front();
        pop_heap(vtx_heap.begin(), vtx_heap.end(), heap_compare());
        vtx_heap.pop_back();
        if(invalid_flag[get<2>(entry)]){
            heapEle result = getApproxMaxGain(triangles[get<2>(entry)], depth);
            max_clique_gains[get<2>(entry)] = make_tuple(result.second, result.first, get<2>(entry));
            vtx_heap.push_back(max_clique_gains[get<2>(entry)]);
            push_heap(vtx_heap.begin(), vtx_heap.end(), heap_compare());
            invalid_flag[get<2>(entry)] = false;
            vtx new_v = result.second;
            size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
            face_store[num_faces-1+max_face_num*new_v] = get<2>(entry);
        }
        else{
            break;
        }
    }
    return entry;
}



