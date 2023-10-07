CC = g++
CCFLAGS = -O3 -Wall  -std=c++17 -mcx16 -mavx2 -flto -march=native #-fopt-info-vec-missed   #-DPROFILE -DPARLAY_SEQUENTIAL# -fcilkplus -DCILKP -DPARLAY_CILK 
CCFLAGS2 = -g -Wall  -std=c++17 -mcx16 -mavx2 -flto -march=native
PLFLAGS = -pthread $(SDSLLF) -flto #-fcilkplus 
PARLAY_DIR=../parlaylib/include
#BOOST_ROOT="../../../../../Downloads/boost_1_82_0"
HIGHWAY_ROOT="../../highway"
GTEST_ROOT="../../highway/googletest-src/googletest/include"


INCLUDE=-I$(PARLAY_DIR) -I$(BOOST_ROOT) -I$(HIGHWAY_ROOT) -I$(GTEST_ROOT) #-L$(HIGHWAY_ROOT) -lhwy

SRC = tmfg.cpp partmfg.cpp partmfg_double.cpp 
NNOBJECTS = $(SRC:.cpp=.o) 



tmfg: $(NNOBJECTS) Makefile
	$(CC) $(PLFLAGS) $(NNOBJECTS) -o $@ -L$(HIGHWAY_ROOT) -lhwy_contrib -lhwy_test -lhwy -L$(HIGHWAY_ROOT)/lib -lgtest -lgtest_main 

tmfg.o: $(SRC) dbht.h linkage.h chain.h partmfg.h
ifdef highway
	$(CC) $(INCLUDE) $(CCFLAGS) -c $(SRC) -D HIGHWAY_MAKE=1 
else
	$(CC) $(INCLUDE) $(CCFLAGS) -c $(SRC)
endif


dynamic: tmfg_dynamic.o
	$(CC) $(PLFLAGS) tmfg_dynamic.o partmfg_double.o -o faster

tmfg_dynamic.o: 
	$(CC) $(INCLUDE) $(CCFLAGS) -c dynamic/tmfg_dynamic.cpp partmfg_double.cpp 


debug: CCFLAGS= -g -O2 -Wall -std=c++17 -DDEBUG -mcx16 -mavx512f -march=native -fno-inline
# -fcilkplus  -DCILKP -DPARLAY_CILK -DPARLAY_USE_STD_ALLOC  -DPARLAY_SEQUENTIAL 

debug: tmfg

openstack: CCFLAGS = -O3 -Wall  -std=c++17 -mcx16 -mavx2 -flto -march=broadwell

openstack: tmfg

openstack: CCFLAGS = -O3 -Wall  -std=c++17 -mcx16 -mavx2 -flto -march=skylake


ec2: tmfg

clean:
	rm -rf tmfg faster *.o
