#!/bin/bash

# Display particle identities
echo "Particle identities:"
echo "  UPQ  = upbar"
echo "  upq  = up"
echo "  glu  = gluon"
echo "  gh   = Ghost"
echo "  GH   = anti-ghost"
echo ""

# Prompt for process input
read -p "Enter the process (e.g., 'glu glu > upq UPQ'): " process

# Loop to ensure the user enters a unique filename
while true; do
    # Prompt for output filename
    read -p "Enter the output filename (e.g., 'gg2uU_LO.qgraf'): " output_filename

    # Check if the output file exists in the Qgraf directory
    if [ -f "Qgraf/$output_filename" ]; then
        echo "Error: File '$output_filename' already exists in the Qgraf directory. Please try again with a different name."
    else
        break
    fi
done

# Ask the user what kind of process to generate: Born, 1-real, or 1-loop
echo "Choose the type of process:"
echo "  1) Born"
echo "  2) 1-real"
echo "  3) 1-loop"
read -p "Enter your choice (1/2/3): " process_choice

# Parse input particles and output particles
in_particles=$(echo $process | cut -d'>' -f1 | xargs)
out_particles=$(echo $process | cut -d'>' -f2 | xargs)

# Assign momenta to input particles
in_momenta=""
i=1
for particle in $in_particles; do
    in_momenta+="${particle}[p$i], "
    i=$((i + 1))
done

# Assign momenta to output particles starting from the last i
out_momenta=""
for particle in $out_particles; do
    out_momenta+="${particle}[p$i], "
    i=$((i + 1))
done

# Remove the trailing comma and space
in_momenta=${in_momenta%, }
out_momenta=${out_momenta%, }

# Set default values for loops and loop momentum
loops="0"
loop_momentum=""
additional_style=""

# Update qgraf.dat based on the process type
case $process_choice in
    1)
        subdir="born"
        ;;
    2)
        subdir="1real"
        out_particles+=" glu[p$i]"
        out_momenta+=", glu[p$i]"
        subdir="1real"
        ;;
    3)
        subdir="1loop"
        loops="1"
        loop_momentum="loop_momentum=k;"
        additional_style="*style= 'xmldraw.sty';"
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

# Generate qgraf.dat content
qgraf_content=$(cat << EOF
output= '$output_filename';
style= 'form_v1.sty';
$additional_style
model= 'QCDGRHIGGS';
in= $in_momenta;
out= $out_momenta;
loops= $loops;
$loop_momentum
options= notadpole, nosnail, onshell;
*true = vsum[gs, 1, 2];
*false = vsum[kap,2,2];
*false = vsum[ew, 1, 2];
*true = vsum[ew, 1, 2];
*true = vsum[kap, 1, 2];
*true = vsum[kq, 1, 2];

* To remove qcd vertex only at oneloop
*true = vsum[gs, 1, 1];
* To keep qcd vertex only at oneloop
*style= 'xmldraw.sty';
EOF
)

# Write to qgraf.dat
echo "$qgraf_content" > qgraf.dat
echo "qgraf.dat file has been updated."

# Inpput file has been updated now we can execute the code.
cd Qgraf
gfortran qgraf-3.1.3.f && ./a.out
cd ../


# Create a directory based on the first letter of each particle and '2' for '>'
process_dir=$(echo $process | sed 's/\bglu\b/g/g' | sed 's/\bupq\b/u/g' | sed 's/\bUPQ\b/U/g' | sed 's/>/2/g' | tr -d ' ')

# Create main directory if it doesn't exist
if [ ! -d "$process_dir" ]; then
    mkdir "$process_dir"
    echo "Main directory '$process_dir' created successfully."
fi

# Create subdirectory based on the process type
if [ ! -d "$process_dir/$subdir" ]; then
    mkdir "$process_dir/$subdir"
    echo "Subdirectory '$process_dir/$subdir' created successfully."
else
    echo "Subdirectory '$process_dir/$subdir' already exists."
fi

