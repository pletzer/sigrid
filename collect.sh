#!/bin/bash -e

## Usage: collect.sh <NAME>

## program to run (relative to build directory)
prog=tests/testConserveInterp2DBig_CXX
## arguments to pass to the program
args=
## command to run gprof2dot.py
#gprof2dot=gprof2dot.py
gprof2dot=/projects/nesi99999/ANTS-regridding/local/bin/gprof2dot.py
## directory for doing the build (can already exist)
build_dir=build

## check argument
if (( $# < 1 ))
then
    echo "Usage: $0 <UNIQUE_NAME> <ARGS>"
    exit 1
fi
if [ -d "$build_dir/$1.dir" ]
then
    echo "Enter unique name"
    exit 2
fi
if (( $# > 1 ))
then
    args="$2"
fi

## load modules
ml purge
ml CMake foss/2015a VTune

## build the code
mkdir -p build
cd build
CC=gcc CXX=g++ cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO ..
make

## run the program to collect data
amplxe-cl -collect hotspots -result-dir $1.dir -search-dir=$(pwd)/cpp -search-dir=$(pwd)/tests -source-search-dir=$(pwd)/../cpp -source-search-dir=$(pwd)/../tests -- $prog $args

## generate gprof-like report
amplxe-cl -report gprof-cc -result-dir $1.dir -format text -report-output $1.txt

## create the call graph
$gprof2dot -f axe $1.txt | dot -Tpng -o $1.png
