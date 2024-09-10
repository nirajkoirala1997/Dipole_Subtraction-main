#!/bin/bash

# Run FORM script
#tform -w5 mat_amp.frm
form mat_amp.frm

# Execute any other necessary shell scripts
bash removenn.sh

# Create a temporary Mathematica script
echo "<< out.m" > mathbro.m
echo "mat // Simplify" >> mathbro.m

# Run the Mathematica kernel interactively and execute the script
math << EOF
<< mathbro.m
EOF
# Clean up the temporary Mathematica script
rm -f mathbro.m

