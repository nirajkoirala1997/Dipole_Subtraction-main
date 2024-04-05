cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/dipole_standalone
make
./runDipole | tee ../output6.dip
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1/summary/compare
gfortran tee_the_data.f -o tee_the_data.o  
./tee_the_data.o 'output6.dip' 'dipole' 'run_06'
rm -f tee_the_data.o
