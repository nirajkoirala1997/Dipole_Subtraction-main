#!/bin/bash
timestamp=$(date +"%Y%m%d%H%M%S")

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_LO.o"

cd $home_path
cd ../LO_check_DY
make clean && make
./runLO | tee "../trash/broken/output_${timestamp}.LO"


cd $home_path
cd ../summary/compare
./tee_the_data_${timestamp}_LO.o "output_${timestamp}.LO" 'LO'

rm -f "tee_the_data_${timestamp}_LO.o"
