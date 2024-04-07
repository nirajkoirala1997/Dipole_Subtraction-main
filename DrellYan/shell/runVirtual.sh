#!/bin/bash
timestamp=$(date +"%Y%m%d%H%M%S")

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_virtual.o"

cd $home_path
cd ../virtual_standalone
make clean && make
./runVir | tee "../trash/broken/output_${timestamp}.virtual"


cd $home_path
cd ../summary/compare
./tee_the_data_${timestamp}_virtual.o "output_${timestamp}.virtual" 'virtual'

rm -f "tee_the_data_${timestamp}_virtual.o"
