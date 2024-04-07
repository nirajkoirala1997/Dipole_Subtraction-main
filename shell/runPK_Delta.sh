#!/bin/bash

timestamp=$(date +"%Y%m%d%H%M%S")

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}_Delta.o"

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/PK_Delta
make clean && make
./runPK | tee "../trash/broken/output_${timestamp}.PK.Delta"


cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
./tee_the_data_${timestamp}_Delta.o "output_${timestamp}.PK.Delta" 'PK_Delta'

rm -f "tee_the_data_${timestamp}_Delta.o"
