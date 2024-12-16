# Compiler
CC = /usr/bin/clang++

# Platform detection
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Common Compiler Flags
COMMON_FLAGS = -O3 -Wall -std=c++17 -stdlib=libc++ -flto
DEBUG_FLAGS = -g -O2 -Wall -std=c++17 -DDEBUG -fno-inline

# Include paths for Parlay and other dependencies
PARLAY_DIR = /usr/local/include
HIGHWAY_DIR = /usr/local/include
GTEST_DIR = /usr/local/include

# Library paths
LIB_PATH = -L/usr/local/lib

# Architecture-specific flags
ifeq ($(UNAME_M), x86_64)
    # x86 (AVX support)
    ARCH_FLAGS = -mcx16 -mavx2 -march=native
else ifeq ($(UNAME_M), arm64)
    # ARM (Apple Silicon)
    ARCH_FLAGS = -march=armv8.3-a
else
    # Generic fallback
    ARCH_FLAGS = -march=native
endif

# Include paths
INCLUDE = -I$(PARLAY_DIR)/parlay -I$(HIGHWAY_DIR) -I$(GTEST_DIR)

# Parallelism flags
PLFLAGS = -pthread

# Source files
SRC = tmfg.cpp partmfg.cpp partmfg_double.cpp
NNOBJECTS = $(SRC:.cpp=.o)

# Target: tmfg
tmfg: $(NNOBJECTS) Makefile
	$(CC) $(PLFLAGS) $(NNOBJECTS) -o $@ $(LIB_PATH) -lhwy_contrib -lhwy_test -lhwy -lgtest -lgtest_main

# Compile object files
%.o: %.cpp
	$(CC) $(INCLUDE) $(COMMON_FLAGS) $(ARCH_FLAGS) -c $< -o $@

# Debug build
debug: CCFLAGS = $(DEBUG_FLAGS) $(ARCH_FLAGS)
debug: tmfg

# Clean target
clean:
	rm -rf tmfg *.o
