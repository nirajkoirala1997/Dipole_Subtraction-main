#!/bin/bash
timestamp=$(date +"%Y%m%d%H%M%S")
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_dipole.o"

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/dipole_standalone/
make clean && make
./runDipole | tee "../trash/broken/output_${timestamp}.dipole"


cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
./tee_the_data_${timestamp}_dipole.o "output_${timestamp}.dipole" 'dipole'

rm -f "tee_the_data_${timestamp}_dipole.o"
