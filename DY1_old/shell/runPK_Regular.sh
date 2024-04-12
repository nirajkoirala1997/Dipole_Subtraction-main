#!/bin/bash

timestamp=$(date +"%Y%m%d%H%M%S")

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_Regular.o"

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/PK_Regular
make clean && make
./runPK | tee "../trash/broken/output_${timestamp}.PK.Regular"


cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
./tee_the_data_${timestamp}_Regular.o "output_${timestamp}.PK.Regular" 'PK_Regular'

rm -f "tee_the_data_${timestamp}_Regular.o"
