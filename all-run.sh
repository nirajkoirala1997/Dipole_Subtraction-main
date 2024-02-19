#!/bin/bash

# Get the current date and time and store it in a variable
current_datetime=$(date +"%Y-%m-%d_%H-%M")

# Create a file with the current date and time as its name
output_file="summary/$current_datetime.txt"

# Display the current date and time on the terminal
echo "Starting processes at $current_datetime"

# Execute the commands and display the output in the terminal and save it to the file
script -c "
{
    touch "$output_file"
    cd PK_standalone
    make
    ./runPK
    cd ..

    cd virtual_standalone
    make
    ./runVir
    cd ..

    cd dipole_standalone
    make
    ./run.x
    cd ..
}" "$output_file"
sed -i 's/\r$//' "$output_file"
sed -i 's/\x1B\[[0-9;]*m//g' "$output_file"

# Open the output file for further editing
vi "$output_file"

