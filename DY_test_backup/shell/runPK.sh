#!/bin/bash

timestamp=$(date +"%Y%m%d%H%M%S")

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY_Q/summary/compare
gfortran tee_the_data.f -o "tee_the_data_${timestamp}.o"

cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY_Q/PK_standalone/
make clean && make
./runPK | tee "../trash/broken/output_${timestamp}.PK"


cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY_Q/summary/compare
./tee_the_data_${timestamp}.o "output_${timestamp}.PK" 'PK'

rm -f "tee_the_data_${timestamp}.o"
