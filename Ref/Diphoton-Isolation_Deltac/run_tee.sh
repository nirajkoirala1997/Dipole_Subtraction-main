#!/bin/bash


# Begin with this .sh file location
home_path=$(dirname "$0")
cd $home_path

# home path assigned for future reference
home_path=$(pwd)

sed -i '2s/0/1/' nfile.dat


# Compile the file to update summary data
cd summary/compare
gfortran tee_the_data.f -o "tee_the_data.o"

#Everything is ready now start the Executable
cd $home_path
make && ./run.x | tee "summary/output.dat"


# Results are now ready we can combine the data of tee to the output files.
cd $home_path
cd summary/compare
./tee_the_data.o 
rm -f "tee_the_data.o"

cd $home_path
sed -i '2s/1/0/' nfile.dat

