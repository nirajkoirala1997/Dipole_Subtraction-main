#!/bin/bash
cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1

# Define the commands to run in parallel, and add a unique identifier to each command
command1="(./shell/runPK.sh)"
command2="(./shell/runDipole.sh)"
command3="(./shell/runVirtual.sh)"

# Run the commands in parallel using GNU Parallel and capture the output
parallel --line-buffer ::: "$command1" "$command2" "$command3" 
