#!/bin/bash
path=$(pwd)
#cd /home/niraj/OneDrive/WorkSpaceIITG/Generalfiles/integrators/DY1
cd $path
cd ../
# Define the commands to run in parallel, and add a unique identifier to each command
command1="(./shell/runPK_Plus.sh)"
command2="(./shell/runPK_Regular.sh)"
command3="(./shell/runPK_Delta.sh)"
command4="(./shell/runDipole.sh)"
command5="(./shell/runVirtual.sh)"

# Run the commands in parallel using GNU Parallel and capture the output
parallel --line-buffer ::: "$command1" "$command2" "$command3" "$command4" "$command5" 
