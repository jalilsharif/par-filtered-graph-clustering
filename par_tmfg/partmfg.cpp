#include "partmfg.h"
#include "partmfg_add.cpp"
#include "partmfg_update.cpp"
#include "dbht.h"



template<class T, class PROF>
void ParTMFG<T, PROF>::heapifyFace(face i){
    pf->incHeapifyNum();
    pf->incHeapifySize(vertex_num);
    triT t = triangles[i];
    auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
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
        #ifdef HIGHWAY_MAKE
        hwy::VQSort(&z[0], vertex_num, hwy::SortAscending());
        #else
        std::sort(&z[0], &z[0]+vertex_num);
        #endif
        for(vtx v = 0; v < vertex_num; v++){
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
}

// sorts the vertices by correlation with a given vertex
template<class T, class PROF>
void ParTMFG<T, PROF>::heapifyVtx(vtx i){
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
        #ifdef HIGHWAY_MAKE
        hwy::VQSort(&z[0], n, hwy::SortAscending());
        #else
        std::sort(&z[0], &z[0]+n);
        #endif
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

template<class T, class PROF>
void ParTMFG<T, PROF>::removeOneV(vtx v){
    vertex_flag[v] = false;
    if(use_corrs){
        vertex_num = n - peo_ind;
        return;
    }
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

template<class T, class PROF>
void ParTMFG<T, PROF>::filterVbyFlag(){
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

template<class T, class PROF> 
void ParTMFG<T, PROF>::initGainArray(){
        if(!hasUninsertedV()) return;
        initGainArrayHeap();
    }

template<class T, class PROF> 
parlay::sequence<size_t> ParTMFG<T, PROF>::getAllBestVertices(size_t THRESHOLD){
        // get the best triangle for each conflict vertex
        parlay::sort_inplace(make_slice(max_clique_gains).cut(0, THRESHOLD));

        auto tmp_lst = parlay::delayed_seq<size_t>(THRESHOLD, [&](size_t i){ return i; });
        auto insert_list = parlay::filter(tmp_lst, [&](auto i){ 
        return i==THRESHOLD-1 || (get<0>(max_clique_gains[i])!= get<0>(max_clique_gains[i+1])); 
        });
        return insert_list;

    }

template<class T, class PROF> 
parlay::sequence<size_t> ParTMFG<T, PROF>::getBestVertices(size_t THRESHOLD){
        if(THRESHOLD < triangles_ind){
            parlay::sort_inplace(make_slice(max_clique_gains).cut(0, triangles_ind), 
            [&](const gainT &i, const gainT &j){
                if(get<1>(i) == get<1>(j)) return i < j;
                return get<1>(i) > get<1>(j); // the largest at the front after sorting
            });
        }
        return getAllBestVertices(THRESHOLD);

    }

template<class T, class PROF> 
void ParTMFG<T, PROF>::updateGainArray(sequence<size_t> &insert_list){
    if(!hasUninsertedV()) return;
    updateGainArrayHeap(insert_list);
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::updateGainArrayHeap(sequence<size_t> &insert_list){
        // auto in = make_slice(vertex_list).cut(vertex_start, vertex_start+vertex_num);
        auto v_list = sequence<pair<vtx,face>>::uninitialized(insert_list.size());
        parlay::parallel_for(0, insert_list.size(), [&](size_t i) {
            auto entry = max_clique_gains[insert_list[i]];
            v_list[i]= make_pair(get<0>(entry), get<2>(entry));
        });
// optimize this
        parlay::sort_inplace(make_slice(max_clique_gains).cut(0, triangles_ind-2*insert_list.size()), 
            [&](const gainT &i, const gainT &j){
                return get<2>(i) < get<2>(j); // sort by face indices to get the corresponding orders again, sort only the before the insert part
            });
        if(use_corrs){
            // get vertex indices to update
            auto indices = sequence<size_t>::uninitialized(v_list.size());
            for(unsigned long int i = 0; i < v_list.size(); i++){
                indices[i] = vtx_to_face_inds[v_list[i].first];
            }
            // remove duplicate vertices
            size_t max_idx = *parlay::max_element(indices, std::less<size_t>{});
            auto update_set = parlay::hashtable(3 * max_idx * v_list.size(), parlay::hash_numeric<T>());
            parlay::parallel_for(0, (size_t)v_list.size(), [&](size_t j) {
                vtx v = v_list[j].first;
                for(size_t ii = 0; ii < vtx_to_face_inds[v]; ii++) {
                    face i = face_store[ii+max_face_num*v];
                    triT new_tri = triangles[i];
                    vtx v1, v2, v3;
                    tie(v1,v2,v3) = new_tri;
                    update_set.insert(v1);
                    update_set.insert(v2);
                    update_set.insert(v3);

                }
            });
            // update vertices
            sequence<T> vertices = update_set.entries();
            parlay::parallel_for(0, vertices.size(), [&](size_t vt) {
                getMaxCorr(vertices[vt]);
                corr_sorted_list_pointer[vertices[vt]]--;
            });
            // update faces
            parlay::parallel_for(0, (size_t)v_list.size(), [&](size_t j) {
                vtx v = v_list[j].first;
                for(size_t ii = 0; ii < vtx_to_face_inds[v]; ii++) {
                    face i = face_store[ii+max_face_num*v];
                    triT new_tri = triangles[i];
                    heapEle result = getFastMaxGain(new_tri, 1);
                    vtx new_v = result.second;
                    max_clique_gains[i] = make_tuple(new_v, result.first, i);
                    max_gains_array[i] = (float)(result.first);
                    size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
                    face_store[num_faces-1+max_face_num*new_v] = i;
                }
            });
        }
        else{
            size_t triangles_ind_prev = triangles_ind - 2*insert_list.size();
            parlay::parallel_for(0, v_list.size(), [&](size_t j) {
                vtx old_v = v_list[j].first; face tri = v_list[j].second;
                parlay::parallel_for(0, vtx_to_face_inds[old_v], [&](size_t ii) {
                    face i = face_store[ii+max_face_num*old_v];
                    if(i == tri || i >= triangles_ind_prev ){heapifyFace(i);}
                    heapEle result = getMinValidHeapEle(i);
                    vtx new_v = result.second;
                    max_clique_gains[i] = make_tuple(new_v, result.first, i);
                    size_t num_faces = pbbs::write_add(&vtx_to_face_inds[new_v], 1);
                    face_store[num_faces-1+max_face_num*new_v] = i; 
                });
            });
        }

    }

template<class T, class PROF> 
void ParTMFG<T, PROF>::insertMultiple(sequence<size_t> &insert_list, DBHTTMFG<T, PROF> *clusterer){
        // for each vertex, triangle pair in the insert_list, insert the vertex 
        // and update triangle list, vertex_flag, and P
        parlay::parallel_for(0, insert_list.size(), [&](size_t i) {
            auto entry = max_clique_gains[insert_list[i]];
            vtx v = get<0>(entry);
            face tri = get<2>(entry);
            vtx t1,t2,t3;
            tie(t1,t2,t3) = triangles[tri];

            peo[peo_ind + i] = v;
            cliques[peo_ind + i-3] = cliqueT(t1,t2,t3,v);
            insertToP(v, t1, P_ind + i*3);
            insertToP(v, t2, P_ind + i*3 + 1);
            insertToP(v, t3, P_ind + i*3 + 2);
            vertex_flag[v] = false;
            triangles[tri] = triT(t1,t2,v);
            triangles[triangles_ind + i*2] = triT(t2,t3,v);
            triangles[triangles_ind + i*2 + 1] = triT(t1,t3,v);

            size_t num_faces = vtx_to_face_inds[v];
            face_store[num_faces+max_face_num*v] = triangles_ind + i*2;
            face_store[num_faces+1+max_face_num*v] = triangles_ind + i*2 + 1;
            vtx_to_face_inds[v]+=2;
            if(clusterer != nullptr){
                clusterer->updateDegrees(t1,t2,t3,v);
                clusterer->insertOne(tri, triangles_ind + i*2, triangles_ind + i*2 + 1, peo_ind+i-3, v);
            }
        });
        peo_ind += insert_list.size();
        triangles_ind += 2*insert_list.size();
        P_ind += 3*insert_list.size();
        filterVbyFlag();
    }




// records graph info, updates triangles and face-vertex pairs
template<class T, class PROF> 
void ParTMFG<T, PROF>::updateInternal(vtx v, face tri, DBHTTMFG<T, PROF> *clusterer){
    vtx t1,t2,t3;
    tie(t1,t2,t3) = triangles[tri];
    peo[peo_ind] = v;
    cliques[peo_ind-3] = cliqueT(t1,t2,t3,v);
    insertToP(v, t1, P_ind);
    insertToP(v, t2, P_ind + 1);
    insertToP(v, t3, P_ind + 2);
    triangles[tri] = triT(t1,t2,v);
    triangles[triangles_ind] = triT(t2,t3,v);
    triangles[triangles_ind + 1] = triT(t1,t3,v);
    size_t num_faces = vtx_to_face_inds[v];
    face_store[num_faces+max_face_num*v] = triangles_ind;
    face_store[num_faces+1+max_face_num*v] = triangles_ind + 1;


    if(clusterer != nullptr){
        clusterer->updateDegrees(t1,t2,t3,v);
        clusterer->insertOne(tri, triangles_ind, triangles_ind + 1, peo_ind-3, v);
    }

    peo_ind += 1;
    triangles_ind += 2;
    P_ind += 3;
    vtx_to_face_inds[v]+=2;
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::insertOne(DBHTTMFG<T, PROF> *clusterer){ // = nullptr
    // get the best vertex among all candidates
    gainT entry = use_max_gains_heap ? getMaxVtxHeap() : getMaxVtxArray();

    vtx v = get<0>(entry);
    face tri = get<2>(entry);

    updateInternal(v, tri, clusterer);
    removeOneV(v); 
    if(vertex_num == 0){
        return;
    }
    
    if(use_max_gains_heap){
        updateVerticesHeap(v, tri);
    }
    else if(use_corrs){
        updateVerticesCorrs(v, tri);
    }
    else{
        updateVerticesOld(v, tri);
    }
}

template<class T, class PROF> 
void ParTMFG<T, PROF>::init(){
        P = sequence<tuple<vtx,vtx,T>>::uninitialized(3*n);
        max_clique_gains = sequence<gainT>::uninitialized(3*n - 6); // ( vertex,gain, face) tuple<vtx, T, face>;
        max_face_num = 3*n - 6;
        cliques = sequence<cliqueT>::uninitialized(n);
        triangles = sequence<triT>::uninitialized(3*n);
        vertex_list = sequence<vtx>::uninitialized(2*n);
        peo = sequence<vtx>::uninitialized(n);

        vertex_flag = parlay::sequence<bool>(n, true);
        face_store = sequence<size_t>::uninitialized(max_face_num*n);
        vtx_to_face_inds = sequence<size_t>(n, 0);

        invalid_flag = sequence<bool>(3*n, false);

        max_gains_array = sequence<float>::uninitialized(3*n);

        cliques[0] = maxClique();

        vtx t1, t2, t3,t4;
        tie(t1,t2,t3,t4) = cliques[0];
        peo[0] = t1; peo[1] = t2; peo[2] = t3; peo[3] = t4;
        for(int i=0; i<4; ++i){
            vertex_flag[peo[i]] = false;
            for(int j=i+1; j<4; ++j){
                insertToP(peo[i], peo[j], P_ind);
                P_ind++;
            }
        }
        triangles[0] = triT(t1, t2, t3);
        triangles[1] = triT(t1, t2, t4);
        triangles[2] = triT(t1, t3, t4);
        triangles[3] = triT(t2, t3, t4);

        auto all_vertices = parlay::delayed_seq<vtx>(n, [&](size_t i){ return i; });
        vertex_num = parlay::filter_into(all_vertices, vertex_list, [&](vtx i){ 
            return i != t1 && i != t2 && i != t3 && i != t4; 
            });
        vertex_start = 0;
        triangles_ind = 4;
        peo_ind = 4;
        heap_buffer = sequence<heapEle>::uninitialized(3*n*n); 
        corr_heap_buffer = sequence<heapEle>::uninitialized(n*n); 
        corr_list = sequence<vtx>::uninitialized(n*n); 


        if(use_corrs){
            corr_sorted_list_pointer=sequence<size_t>(n, 0);
            sorted_list_pointer=sequence<size_t>(3*n, 0);

            parlay::parallel_for(0, n, [&](size_t v){
                heapifyVtx(v);
            });
        }
        sorted_list_pointer=sequence<size_t>(3*n, 0);

    }


template<class T, class PROF> 
void ParTMFG<T, PROF>::initGainArrayHeap(){
    //vtx_heap = new vector<gainT>;
    vtx_heap = vector<gainT>();
    for(face i = 0; i < triangles_ind; i++) {
        heapifyFace(i);
        triT t = triangles[i];
        // issue: heapifyFace (only necessary when use_max_gains_heap is false) does not work properly with highway
        heapEle result = getMinValidHeapEle(i);
        vtx v = result.second;
        max_clique_gains[i] = make_tuple(v, result.first, i);
        max_gains_array[i] = (float)(result.first);
        size_t num_faces = pbbs::write_add(&vtx_to_face_inds[v], 1);
        face_store[num_faces-1+max_face_num*v] = i;
        if(use_max_gains_heap){
            vtx_heap.push_back(max_clique_gains[i]);
        }
    }
    if(use_max_gains_heap){
        make_heap(vtx_heap.begin(), vtx_heap.end(), heap_compare());
    }
}

template<class T, class PROF> 
cliqueT ParTMFG<T, PROF>::maxClique(){
    auto Wseq = parlay::delayed_seq<T>(n*n, [&](size_t i) {
        return getW(i/n, i%n);
    });
    double W_mean = parlay::reduce(Wseq) / n / n;
    auto W_larger = parlay::delayed_seq<T>(n*n, [&](long i){ return Wseq[i] > W_mean? Wseq[i]:0; });


    sequence<pair<T, vtx>> v_weight = sequence<pair<T, vtx>>(n);
    parlay::parallel_for(0, n, [&](size_t i) {
        v_weight[i] = make_pair(parlay::reduce(make_slice(W_larger).cut(i*n, (i+1)*n)), i);
    });

    parlay::sort_inplace(v_weight); // can improve to 4 max operations

    return cliqueT(get<1>(v_weight[n-1]), get<1>(v_weight[n-2]), get<1>(v_weight[n-3]), get<1>(v_weight[n-4]));
}