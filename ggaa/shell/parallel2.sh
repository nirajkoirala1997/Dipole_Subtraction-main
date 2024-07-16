#!/bin/bash

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../

# Define the commands to run in parallel, and add a unique identifier to each command

commands=(
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
    "(./shell/runPK_Plus.sh)"
)

# Define the new values for the "max # of distribution increment step_size from xq"
distribution_steps=(100 400 700 1000 1300 1600 1900 2200 2500)

# Function to modify the input files
modify_input_files() {
    local index=$1
    local step_size=${distribution_steps[$index]}
	echo $step_size

    # Modify run.machine.dat
    sed -i "5s/.*/${step_size}/" run.machine.dat

    # Modify output_files.dat
    sed -i "s/PK_Plus[0-9]*.dat/PK_Plus${index}.dat/" output_files.dat
}

# Create a temporary directory to hold copies of the modified input files
tmp_dir=$(mktemp -d)
cp run.machine.dat output_files.dat $tmp_dir


# Export the modify_input_files function for parallel
export -f modify_input_files



parallel --line-buffer ::: \
    "(sleep 0; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 0)" \
    "(sleep 1; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 1)" \

## Run the commands in parallel using GNU Parallel
#parallel --line-buffer ::: \
#    "(sleep 0; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 0)" \
#    "(sleep 1; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 1)" \
#    "(sleep 2; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 2)" \
#    "(sleep 3; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 3)" \
#    "(sleep 4; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 4)" \
#    "(sleep 5; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 5)" \
#    "(sleep 6; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 6)" \
#    "(sleep 7; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 7)" \
#    "(sleep 8; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 7)"
# Clean up the temporary directory
#rm -rf $tmp_dir

#    "(sleep 0; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 0; ${commands[0]})" \
#    "(sleep 1; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 1; ${commands[1]})" \
#    "(sleep 2; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 2; ${commands[2]})" \
#    "(sleep 3; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 3; ${commands[3]})" \
#    "(sleep 4; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 4; ${commands[4]})" \
#    "(sleep 5; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 5; ${commands[5]})" \
#    "(sleep 6; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 6; ${commands[6]})" \
#    "(sleep 7; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 7; ${commands[7]})" \
#    "(sleep 8; cp $tmp_dir/run.machine.dat $tmp_dir/output_files.dat .; modify_input_files 7; ${commands[8]})"
