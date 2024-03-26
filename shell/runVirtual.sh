cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/virtual_standalone/
make
./runVir | tee ../output.vir
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o tee_the_data.o  
./tee_the_data.o 'output.vir' 'virtual'
rm -f tee_the_data.o
