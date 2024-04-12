#!/bin/bash

timestamp=$(date +"%Y%m%d%H%M%S")

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_Plus.o"

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/PK_Plus
make clean && make
./runPK | tee "../trash/broken/output_${timestamp}.PK.Plus"


cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
./tee_the_data_${timestamp}_Plus.o "output_${timestamp}.PK.Plus" 'PK_Plus'

rm -f "tee_the_data_${timestamp}_Plus.o"
