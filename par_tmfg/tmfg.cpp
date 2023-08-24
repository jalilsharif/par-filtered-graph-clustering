#include <iostream>

// #include "dbht.h"
#include "partmfg_double.h"
#include "profiler.h"
#include "IO.h"

#include <fstream>
#include <string>
#include <getopt.h>

ostream* IO::time_output = &cout;

void runDBHT(SymM<double> *W, SymM<double> *D, size_t n, size_t THRESHOLD, string method, bool use_corrs, bool use_gains_heap, bool use_highway, bool exact_apsp, string dsname = ""){//, bool use_gains_heap = false){
    //ofstream outfile;
    //outfile.open("timedata.txt", std::ios_base::app);

    (*IO::time_output) << "====" << endl;
    (*IO::time_output) << "threshold: " << THRESHOLD << endl;
    (*IO::time_output) << "method: " << method << endl;
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
    (*IO::time_output) << "init total: "<< t.next() << endl;

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
    (*IO::time_output) << "tmfg total: "<< t.next() << endl;
    (*IO::time_output) << "round: " << round << endl;
    computer.computeCost();
    pf.report();
    t.next();
    clusterer->APSP();
    (*IO::time_output) << "APSP total: "<< t.next() << endl;
    clusterer->computeDirection();
    (*IO::time_output) << "direction total: "<< t.next() << endl;
    clusterer->nonDiscreteClustering();
    (*IO::time_output) << "non-discrete total: "<< t.next() << endl;
    clusterer->assignToConvergingBubble(); // need to test
    (*IO::time_output) << "discrete total: "<< t.next() << endl;
    (*IO::time_output) << "num cluster: "<< clusterer->nc << endl;
    clusterer->assignToBubble(); // need to test
    (*IO::time_output) << "bubble total: "<< t.next() << endl;
    clusterer->buildHierarchy();
    (*IO::time_output) << "hierarchy total: "<< t.next() << endl;


    if(method == "exact" || method == "naive"){
        computer.outputP("./outputs/Ps/" + dsname + "-" + method + "-P-1");
        clusterer->outputDendro("./outputs/Zs/" + dsname + "-" + method + "-Z-1" );
    }else{
        computer.outputP("./outputs/Ps/" + dsname + "-" + method + "-P-" + to_string(THRESHOLD) );
        clusterer->outputDendro("./outputs/Zs/" + dsname + "-" + method + "-Z-" + to_string(THRESHOLD));
    }
    (*IO::time_output) << endl;
    cout << "done" << endl;

}

//ostream* time_output = NULL;





int main(int argc, char *argv[]) {


char* filename; // path to correlation data
string dsname; // path to dendrogram output directory
bool dist_file = false; // whether to use distance file or sqrt(2*(1-w)) by default
char* distance_filename; // path to distance file
size_t n; // number of items
string method = "exact"; // method for TMFG construction
size_t THRESHOLD = 0; // number of nodes inserted per iteration if method = "prefix"
int round = 1; // number of times to execute 
bool use_gains_heap = true; // whether to use a heap for TMFG construction, only relevant if use_corrs = true
bool use_corrs = true; // whether to use the new TMFG construction method
bool use_highway = false; // whether to use highway for sorting
bool exact_apsp = false; // whether to use exact or approximate APSP



int opt;
ofstream outfile2;

  while (1) {
        static struct option long_options[] = {
            {"data-file", required_argument, 0, 'f'},
            {"out-file", required_argument, 0, 'o'},
            {"n", required_argument, 0, 'n'},
            {"help", no_argument, 0, 'h'},
            {"dist-file", required_argument, 0, 'd'},
            {"time-output", required_argument, 0, 't'},
            {"prefix", required_argument, 0, 'p'},
            {"rounds", required_argument, 0, 'r'},
            {"no-heap", no_argument, 0, 'N'},
            {"orig-tmfg", no_argument, 0, 'O'},
            {"use-highway", no_argument, 0, 'W'},
            {"exact-apsp", no_argument, 0, 'A'},
            {0, 0, 0, 0}};

        opt = getopt_long(argc, argv, "f:o:n:hd:t:p:r:NOWA", long_options, NULL);
        if (opt == -1)
            break;

        switch (opt) {
        case 'f':
            filename = optarg;
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
        case 'r':
            round = atoi(optarg);
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
        case 'h':
            cout << "usage:\n";
            cout << "Required options, all take arguments:\n";
            cout << "   --data-file, -f:    path to correlation data\n";
            cout << "   --out-file, -o:     path to dendrogram output directory\n";
            cout << "   -n:                 number of items to cluster\n";
            cout << "Other options that take arguments:\n";
            cout << "   --dist-file, -d:    path to input dissimilarity matrix. If none specified, uses sqrt(2(1-s))\n";
            cout << "   --time-output, -t:  path to output file for timing data. If none specified, uses standard output\n";
            cout << "   --prefix, -p:       prefix size to insert in each TMFG iteration, uses exact method (prefix=1) if not specified\n";
            cout << "   --rounds, -r:       number of times to run the algorithm, defaults to 1\n";
            cout << "Other options that do not take arguments:\n";
            cout << "   --no-heap, -N:      scans the array of vertices each iteration during TMFG construction\n";
            cout << "   --orig-tmfg, -O:    uses the original TMFG construction method\n";
            cout << "   --use-highway, -W:  uses highway instead of std::sort and parlay sort\n";
            cout << "   --exact-apsp, -A:   uses exact all-pairs shortest paths instead of approximate\n";
            cout << "   --help, -h:         prints help message; no other options needed\n";
            return(0);
            break;
        case '?':
            break;
        default:
            abort();
        }
  }

#ifndef HIGHWAY_MAKE
if(use_highway){
    cout << "Warning: --use-highway option selected but highway is disabled.\n";
}
use_highway = false;
#endif //HIGHWAY_MAKE

(*IO::time_output) << dsname << ' ';

(*IO::time_output) << use_corrs << use_gains_heap<<use_highway<<endl;
(*IO::time_output) << "workers: " << parlay::num_workers() << endl;
timer t2;t2.start();
SymM<double> W = IO::readSymMatrixFromFile<double>(filename, n); 
W.setDiag(1);
for(int r=0;r<round;++r){
    if(dist_file){
        SymM<double> D = IO::readSymMatrixFromFile<double>(distance_filename, n);
        (*IO::time_output) << "read: " << t2.next() << endl;
        runDBHT(&W, &D, n, THRESHOLD, method, use_corrs, use_gains_heap, use_highway, exact_apsp, dsname);
        D.free_matrix();
    }else{
        (*IO::time_output) << "read: " << t2.next() << endl;
        runDBHT(&W, nullptr, n, THRESHOLD, method, use_corrs, use_gains_heap, use_highway, exact_apsp, dsname);
    }
}

}
