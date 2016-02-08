CC=mpic++.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=$(HOME)/workspace/sgpp-trunk

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ -I$(SGPP_DIR)/combigrid/src/

LDIR=$(LD_SGPP)
INC=$(INC_SGPP)

LIBS=-lsgppdistributedcombigrid -lsgppcombigrid -lsgppbase -lboost_serialization

all: combi_example 

combi_example: combi_example.cpp TaskExample.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_example combi_example.cpp $(LIBS)

clean:
	rm -f *.o out/* combi_example
