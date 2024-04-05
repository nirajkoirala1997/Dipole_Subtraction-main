cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/virtual_standalone/
make
./runVir | tee ../output6.vir
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o tee_the_data.o  
./tee_the_data.o 'output6.vir' 'virtual' 'run_06'
rm -f tee_the_data.o
