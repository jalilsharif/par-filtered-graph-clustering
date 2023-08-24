#include "dbht.h"

using namespace std;

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::outputDirections(){//string filename, size_t _n=0){
        for(size_t i=0;i<n-3;i++){
            cout << directions[i] << endl;
        }
    } 

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::outputDefaultClustering(string filename, size_t _n){
        if(_n==0)_n = n;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << flatClustering[i] << endl;
        }
        file_obj.close();
    } 

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::outputBubbleAssign(string filename, size_t _n){
        if(_n==0)_n = n;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            file_obj << bbMember[i] << endl;
        }
        file_obj.close();
    } 

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::outputDendro(string filename, size_t _n){
        if(_n==0)_n = n-1;
        ofstream file_obj;
        file_obj.open(filename); 
        for(size_t i=0;i<_n;i++){
            dendro[i].print(file_obj);
        }
        file_obj.close();
    } 

template<class T, class PROF> 
void DBHTTMFG<T, PROF>::outputLabels(string filename, int time){
        ofstream file_obj;
        file_obj.open(filename, std::ofstream::out|std::ios_base::app);
        if (file_obj.fail()) { 
            cout<<"failed...";
        }
        file_obj << time << endl;
        for(size_t i=0;i<n;i++){
            file_obj << cluster_labels[i] << ' ';
        }
        file_obj << endl;
        file_obj.close();
    } 

