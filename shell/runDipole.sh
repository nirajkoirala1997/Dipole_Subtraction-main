cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/dipole_standalone
make
./runDipole | tee ../output.dip
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o tee_the_data.o  
./tee_the_data.o 'output.dip' 'dipole'
rm -f tee_the_data.o
