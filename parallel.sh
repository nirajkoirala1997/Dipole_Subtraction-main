#!/bin/bash

# Modify the Fortran source file
sed -i 's/nprn[[:space:]]*=[[:space:]]*1/nprn = 0/' "vsup.f"

# Define the commands to run in parallel, and add a unique identifier to each command
command1="(cd dipole_standalone/ && make  && ./runDipole)"
command2="(cd PK_standalone/ && make && ./runPK)"
command3="(cd virtual_standalone/ && make && ./runVir)"

# Run the commands in parallel using GNU Parallel and capture the output
parallel --line-buffer ::: "$command1" "$command2" "$command3" | tee output.log

# Revert the modification in the Fortran source file
sed -i 's/nprn[[:space:]]*=[[:space:]]*0/nprn = 1/' "vsup.f"

