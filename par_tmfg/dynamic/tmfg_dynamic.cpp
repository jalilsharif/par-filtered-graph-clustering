#include "tmfg_dynamic.hpp"
#include "dynamic_io.cpp"

#include "../partmfg_double.h"
#include "../profiler.h"
#include "../IO.h"

#include <getopt.h>



ostream* IO::time_output = &cout;

template <class T>
sequence<int> dynamic_TMFG<T>::runDBHT(SymM<double> *W, SymM<double> *D, size_t n, size_t THRESHOLD, string method, bool use_corrs, bool use_gains_heap, bool use_highway, bool exact_apsp, int num_clusters, int time, bool relabel, sequence<int> prev_labels, bool verbose, string dsname){//, bool use_gains_heap = false){
    //ofstream outfile;
    //outfile.open("timedata.txt", std::ios_base::app);

    //(*IO::time_output) << "====" << endl;
    //(*IO::time_output) << "threshold: " << THRESHOLD << endl;
    //(*IO::time_output) << "method: " << method << endl;
    #ifdef PROFILE
        auto pf = Profiler();
    #else
        auto pf = DummyProfiler();
    #endif
    timer t2;t2.start();
    ParTMFGD computer = ParTMFGD(W, n, &pf, use_corrs, use_gains_heap, use_highway); 
    timer t;t.start();
    computer.init();   
    computer.initGainArray();                 pf.setInitTime(t2.next());
    auto clusterer = new ParDBHTTMFGD(computer.cliques.data(), computer.triangles.data(), n, computer.W, computer.P.data(), D, &pf, exact_apsp);
    int round=0;
    if(verbose){
        (*IO::time_output) << "init total: "<< t.next() << endl;
    }

    if(method == "prefix"){ //get best gain by scanning vertex list
        while(computer.hasUninsertedV()){
                size_t round_THRESHOLD = min(THRESHOLD, computer.getTrianglesNum());
                auto insert_list = computer.getBestVertices(round_THRESHOLD);   pf.incVTime(t2.next());
                computer.insertMultiple(insert_list, clusterer);                 pf.incInsertTime(t2.next());
                computer.updateGainArray(insert_list);                          pf.incUpdTime(t2.next());
    #ifdef DEBUG
            computer.checkTriangles();
    #endif
                round++;
        } //while end
    }else if(method == "exact"){ //use exact
        while(computer.hasUninsertedV()){
            computer.insertOne(clusterer);
            round++;
        } //while end
    }else if(method == "naive"){ //naive method
        while(computer.hasUninsertedV()){
            auto insert_list = computer.getAllBestVertices(computer.getTrianglesNum()); pf.incVTime(t2.next());
            computer.insertMultiple(insert_list, clusterer);                             pf.incInsertTime(t2.next());
            computer.initGainArray();                                                   pf.incUpdTime(t2.next());
            round++;
        } //while end
    }
    // compute the total gain here
    if(verbose){
        (*IO::time_output) << "tmfg total: "<< t.next() << endl;
        (*IO::time_output) << "round: " << round << endl;
        computer.computeCost();
    }
    pf.report();
    t.next();
    clusterer->APSP();
    if(verbose){
        (*IO::time_output) << "APSP total: "<< t.next() << endl;
    }
    clusterer->computeDirection();
    if(verbose){
        (*IO::time_output) << "direction total: "<< t.next() << endl;
    }
    clusterer->nonDiscreteClustering();
    if(verbose){
        (*IO::time_output) << "non-discrete total: "<< t.next() << endl;
    }
    clusterer->assignToConvergingBubble(); // need to test
    if(verbose){
        (*IO::time_output) << "discrete total: "<< t.next() << endl;
        (*IO::time_output) << "num cluster: "<< clusterer->nc << endl;
    }
    clusterer->assignToBubble(); // need to test
    if(verbose){
        (*IO::time_output) << "bubble total: "<< t.next() << endl;
    }
    clusterer->buildHierarchy();
    if(verbose){
        (*IO::time_output) << "hierarchy total: "<< t.next() << endl;
    }

    if(verbose){
        if(method == "exact" || method == "naive"){
            computer.outputP("outputs/Ps/" + dsname + "-" + method + "-P-1");
            clusterer->outputDendro("outputs/Zs/" + dsname + "-" + method + "-Z-1"  );
        }else{
            computer.outputP("outputs/Ps/" + dsname + "-" + method + "-P-" + to_string(THRESHOLD) );
            clusterer->outputDendro("outputs/Zs/" + dsname + "-" + method + "-Z-" + to_string(THRESHOLD));
        }
    }


    if(relabel){
        for(int i=0;i<n;i++){
                clusterer->prev_labels[i]=prev_labels[i];
        }
    }
    clusterer->clusterLabels(num_clusters);
    if(relabel){
        clusterer->relabel(num_clusters);
    }
    clusterer->outputLabels("dynamic/dynout.txt",time);
    
    //(*IO::time_output) << endl;
    
    return clusterer->cluster_labels;

}


template <class T>
void dynamic_TMFG<T>::tick(bool start, bool stop, int i, int clusters, bool verbose, string method, int THRESHOLD, bool use_corrs, bool use_gains_heap, bool use_highway, bool exact_apsp, string dsname)
{
    corrs.update_corr_matrix();
    //string dsname = "ECG5000";
    //string method = "exact";
    
    if(true){ // replace with better heuristic
        cout << "Running DBHT at time " << i <<'\n';
        prev_labels=runDBHT(&(corrs.corr_matrix), nullptr, n, 0, method, use_corrs, use_gains_heap, use_highway, exact_apsp, clusters, i, !start, prev_labels, verbose, dsname);
        cout << "done\n";
        recluster_ct++;
    }
       
    if(!stop){
        corrs.advance();
    }
    if(verbose){
        cout<<"tick\n";
    }
}


int main(int argc, char *argv[]) {
    bool twofiles = false;
    string file1 = "../../../../../large_files/UCRArchive_2018/ECG5000/ECG5000_TEST.tsv";
    string file2 = "../../../../../large_files/UCRArchive_2018/ECG5000/ECG5000_TRAIN.tsv";
    string fname = "ecg.tsv";
    int window_size = 10;
    int series_len = 140;
    string labels = "dynamic/dynout.txt";
    int clusters = 5;
    bool verbose = false; // whether to output everything or just labels
    string dsname = "ECG5000"; // path to dendrogram output directory
    size_t n = 5000;
    string method = "exact";
    size_t THRESHOLD = 0;
    int round = 1;
    bool use_gains_heap = true;
    bool use_corrs = true;
    bool use_highway = false;
    bool exact_apsp = false;

    bool dist_file = false;
    string distance_filename = "";



    int opt;
    ofstream outfile2;

    while(1){
        static struct option long_options[]{
            {"two-files", no_argument, 0, '2'},
            {"file1", required_argument, 0, 'q'},
            {"file2", required_argument, 0, 'z'},
            {"data-file", required_argument, 0, 'f'},
            {"window-size", required_argument, 0, 'w'},
            {"length", required_argument, 0, 's'},
            {"clusters", required_argument, 0, 'c'},
            {"verbose", no_argument, 0, 'v'},
            {"label-file", required_argument, 0, 'l'},
            {"info-name", required_argument, 0, 'o'},
        	{"n", required_argument, 0, 'n'},
        	{"help", no_argument, 0, 'h'},
        	{"dist-file", required_argument, 0, 'd'},
        	{"time-output", required_argument, 0, 't'},
        	{"prefix", required_argument, 0, 'p'},
        	{"no-heap", no_argument, 0, 'N'},
        	{"orig-tmfg", no_argument, 0, 'O'},
        	{"use-highway", no_argument, 0, 'W'},
        	{"exact-apsp", no_argument, 0, 'A'},


            {0,0,0,0}
        };

        opt = getopt_long(argc, argv, "2q:z:f:w:s:c:vl:o:n:hd:t:p:NOWA", long_options, NULL);

        if (opt == -1)
        	break;

    	switch (opt) {
    	case 'f':
        	fname = optarg;
        	break;
    	case 'o':
        	dsname = optarg;
        	break;
    	case 'n':
        	n = atoi(optarg);
        	break;
    	case 'd':
        	dist_file = true;
        	distance_filename = optarg;
        	break;
    	case 't':
        	outfile2.open(optarg, std::ios_base::app);
        	IO::time_output = &outfile2;
        	break;
    	case 'p':
        	method = "prefix";
        	THRESHOLD = atoi(optarg);
        	break;
        case 'N':
        	use_gains_heap = false;
        	break;
    	case 'O':
        	use_corrs = false;
        	break;
    	case 'W':
        	use_highway = true;
        	break;
    	case 'A':
        	exact_apsp = true;
        	break;
        case 'v':
            verbose = true;
            break;
        case '2':
            twofiles = true;
            break;
        case 'q':
            file1 = optarg;
            break;
        case 'z':
            file2 = optarg;
            break;
        case 'w':
            window_size = atoi(optarg);
            break;
        case 's':
            series_len = atoi(optarg);
            break;
        case 'l':
            labels = optarg;
            break;
        case 'c':
            clusters = atoi(optarg);
            break;
        case '?':
            break;
        default:
            abort();
        }
        
    }


    std::ofstream clearer;
    clearer.open(labels, std::ofstream::out | std::ofstream::trunc);
    clearer.close();
    //int num_points = 140;
    //int window = 10;

    dynamic_TMFG<double> d = dynamic_TMFG<double>(n, series_len, window_size, twofiles, file1, file2, fname);

    
    for(int time=window_size-1; time<series_len; time++){
        bool start = (time==window_size-1);
        bool stop = (time==series_len-1);
        d.tick(start, stop, time, clusters, verbose, method, THRESHOLD, use_corrs, use_gains_heap, use_highway, exact_apsp, dsname);
    }
}