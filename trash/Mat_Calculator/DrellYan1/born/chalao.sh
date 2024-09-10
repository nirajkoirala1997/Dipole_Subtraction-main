#!/bin/bash

# Run FORM script
tform -w6 mat_amp.frm > out.m

# Execute any other necessary shell scripts
sed -i '1,4d' out.m
sed -i 's/p1.nv/(&)/g' out.m
sed -i 's/p2.nv/(&)/g' out.m
sed -i 's/p3.nv/(&)/g' out.m
sed -i 's/p4.nv/(&)/g' out.m
sed -i 's/p5.nv/(&)/g' out.m
sed -i 's/mat =/mat =(/g' out.m
sed -i 's/;/);/g' out.m

# Create a temporary Mathematica script
echo "<< out.m" > mathbro.m
echo "mat // Simplify " >> mathbro.m
#echo "mat/.{s12->s,s13->-t,s23->-u} // Simplify " >> mathbro.m

# Run the Mathematica kernel interactively and execute the script
math << EOF
<< mathbro.m
EOF
# Clean up the temporary Mathematica script
rm -f mathbro.m

