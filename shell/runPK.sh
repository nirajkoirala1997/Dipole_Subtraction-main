cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/PK_standalone/
make clean && make
./runPK | tee  ../output.PK
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o tee_the_data.o  
./tee_the_data.o 'output.PK' 'PK' 'run_06' 
rm -f tee_the_data.o
