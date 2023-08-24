#include "parlay/primitives.h"
#include "parlay/sequence.h"


#include "linkage_types.h"
#include "dbht.h"

#include <queue>
#include <vector>

using parlay::sequence;

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::clusterLabels(int num_clusters){
    //sequence<dendroLine> dendro_cpy = sequence<dendroLine>::uninitialized(n-1);
    //std::memcpy(dendro_cpy.data(), dendro, sizeof(dendroLine)*(n-1));
    //parlay::internal::quicksort(dendro_cpy, compareHeightsReverse);
    sequence<int> cluster_roots = sequence<int>::uninitialized(num_clusters);
    //https://stackoverflow.com/questions/14902876/indices-of-the-k-largest-elements-in-an-unsorted-length-n-array
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::less<std::pair<double, int>> > top_nodes_queue;
    for (int i = 0; i < n-1; i++){
        top_nodes_queue.push(std::pair<double, int>(dendro[i].height, i+n));
    }
    // what if there are too many clusters compared to nodes?
    sequence<int> top_nodes = sequence<int>::uninitialized(num_clusters);
    for(int i = 0; i < num_clusters; i++){
        top_nodes[i] = top_nodes_queue.top().second;
        top_nodes_queue.pop();
    }
    // cut off tree to get individual clusters
    cluster_roots[0] = top_nodes[0];
    for(int i = 1; i < num_clusters; i++){
        int cur_node = top_nodes[i-1];
        cluster_roots[i] = getDendroNode(cur_node).id2;
        for(int j=0; j<i; j++){
            if(cluster_roots[j]==cur_node){
                cluster_roots[j] = getDendroNode(cur_node).id1;
                break;
            }
        }
    }

    cluster_labels = sequence<int>(n, -1);

    std::function<void(int, int)> process_cluster = [&](int root, int label){
        if(root < n){
            cluster_labels[root] = label;
        }
        else{
            process_cluster(getDendroNode(root).id1, label);
            process_cluster(getDendroNode(root).id2, label);
        }
    };

    
    for(int label=0; label<num_clusters; label++){
        process_cluster(cluster_roots[label], label);
    }


}

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::relabel(int num_clusters){
    sequence<int> relabel_dict = sequence<int>(num_clusters, -1);
    sequence<int> relabel_inverse = sequence<int>(num_clusters, -1);
    sequence<int> pair_counter = sequence<int>(num_clusters*num_clusters, 0);
    sequence<int> prev_single_counter = sequence<int>(num_clusters, 0);
    sequence<int> cur_single_counter = sequence<int>(num_clusters, 0);
    for(int idx=0; idx<n-1; idx++){
        int cur_label = cluster_labels[idx];
        int prev_label = prev_labels[idx];
        cur_single_counter[cur_label]++;
        prev_single_counter[prev_label]++;
        pair_counter[cur_label*num_clusters+prev_label]++;
    }
    for(int cur_label=0; cur_label<num_clusters; cur_label++){
        double max_sim = 0;
        int best_base = 0;
        for(int base=0; base<num_clusters; base++){
            if(relabel_inverse[base] != -1){
                continue;
            }
            // div by 0?
            double new_sim = (double)pair_counter[cur_label*num_clusters+base]/(double)(prev_single_counter[base]+cur_single_counter[cur_label] - pair_counter[cur_label*num_clusters+base]);
            if(new_sim > max_sim){
                best_base = base;
                max_sim = new_sim;
            }
        }
        relabel_dict[cur_label] = best_base;
        relabel_inverse[best_base] = cur_label;
    }
    for(int idx=0; idx<n; idx++){
        cluster_labels[idx] = relabel_dict[cluster_labels[idx]];
    }
}