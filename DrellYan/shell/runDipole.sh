#!/bin/bash
timestamp=$(date +"%Y%m%d%H%M%S")

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_dipole.o"
#
cd $home_path
cd  ../dipole_standalone
make clean && make
./runDipole | tee "../trash/broken/output_${timestamp}.dipole"
#
#
cd $home_path
cd ../summary/compare
./tee_the_data_${timestamp}_dipole.o "output_${timestamp}.dipole" 'dipole'
#
rm -f "tee_the_data_${timestamp}_dipole.o"
