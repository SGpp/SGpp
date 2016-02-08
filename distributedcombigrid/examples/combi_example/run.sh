#!/bin/bash
export LD_LIBRARY_PATH=:/home/sccs/workspace/combi/lib/sgpp:/home/heenemo/workspace/sgpp-trunk/lib/sgpp
mpirun.mpich -n 5 ./combi_example
