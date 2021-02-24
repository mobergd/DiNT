#!/bin/bash

FC=ftn

#rm -r build/*
if [ ! -d "build/" ]; then
  mkdir build/
fi
cd build/
#cmake -DCMAKE_Fortran_FLAGS=" -O3 -g -L/home/ajasper/KTP/dint/sprng/lib -L/usr/lib64 -L/home/moberg/lib -lmlfg" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../ ../
cmake -DCMAKE_Fortran_FLAGS="-mkl" -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_INSTALL_PREFIX=../ ../
#cmake -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_INSTALL_PREFIX=../ ../
make VERBOSE=1
make install

#cd src/
#make MakefileOLD
