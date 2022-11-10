#!/bin/bash
#

export SPGLIB_ROOT='<this_code_path>/external/spglib/install'

# Compile the program
g++ ../src/spglib_2d_tovasp.cpp -o spglib_2d_tovasp.x \
  -I${SPGLIB_ROOT}/include \
  -L${SPGLIB_ROOT}/lib/ \
  -lsymspg

# Execute it!
export LD_LIBRARY_PATH=${SPGLIB_ROOT}/lib/:${LD_LIBRARY_PATH}
./spglib_2d_tovasp.x -f POSCAR -p 1E-7 $*
