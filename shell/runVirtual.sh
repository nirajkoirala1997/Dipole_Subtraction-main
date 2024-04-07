#!/bin/bash
timestamp=$(date +"%Y%m%d%H%M%S")
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_virtual.o"

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/virtual_standalone/
make clean && make
./runVir | tee "../trash/broken/output_${timestamp}.virtual"


cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
./tee_the_data_${timestamp}_virtual.o "output_${timestamp}.virtual" 'virtual'

rm -f "tee_the_data_${timestamp}_virtual.o"
