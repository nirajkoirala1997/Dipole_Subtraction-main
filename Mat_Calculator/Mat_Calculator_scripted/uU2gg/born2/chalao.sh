#!/bin/bash

# Run FORM script
#tform -w5 mat_amp.frm
./rum.sh

# Execute any other necessary shell scripts
bash removenn.sh

# Create a temporary Mathematica script
echo "<< out2.m" > mathbro.m
echo "mat /.u->-s-t // Simplify" >> mathbro.m

# Run the Mathematica kernel interactively and execute the script
math << EOF
<< mathbro.m
EOF
# Clean up the temporary Mathematica script
rm -f mathbro.m

