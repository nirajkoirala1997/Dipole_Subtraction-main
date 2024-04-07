#!/bin/bash
timestamp=$(date +"%Y%m%d%H%M%S")

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_PK.Delta.o"

cd $home_path
cd ../PK_Delta
make clean && make
./runPK | tee "../trash/broken/output_${timestamp}.PK.Delta"


cd $home_path
cd ../summary/compare
./tee_the_data_${timestamp}_PK.Delta.o "output_${timestamp}.PK.Delta" 'PK_Delta'

rm -f "tee_the_data_${timestamp}_PK.Delta.o"
