cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/PK_test3/
make
./runPK | tee  ../output.PK
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o tee_the_data.o  
./tee_the_data.o 'output.PK' 'PK' 
rm -f tee_the_data.o
