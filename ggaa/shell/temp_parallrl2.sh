#!/bin/bash

home_path=$(dirname "$0")
cd $home_path
home_path=$(pwd)
cd ../

output_dir=$(sed -n '7p' run.machine.dat | awk '{print $1}')

output_summary_file="summary/${output_dir}/Plus_all.dat"
temp_file="temp_summary.dat"

# Ensure the summary file is empty before we start
> $output_summary_file
> $temp_file

# Define the Q values
distribution_steps=(100 400 700 1000 1300 1600 1900 2200 2500)

# Iterate over the distribution steps and process the corresponding output files
for i in ${!distribution_steps[@]}; do
    q_value=${distribution_steps[$i]}
    file="summary/${output_dir}/PK_Plus${i}.dat"
    
    # Check if the file exists
    if [[ -f $file ]]; then
        # Read the first line and extract the required values
        first_line=$(head -n 1 $file)
        # Add the first line's data to the temporary summary file
        echo $first_line | awk -v q_value=$q_value '{printf "%10d %25.15e %25.15e\n", q_value, $2, $3}' >> $temp_file
        # Append the rest of the file's content to the output summary file
        tail -n +2 $file >> $output_summary_file
    else
        echo "Warning: File $file not found."
    fi
done

# Prepend the processed first lines to the output summary file
cat $temp_file $output_summary_file > "temp_all.dat" && mv "temp_all.dat" $output_summary_file

# Remove the temporary file
rm -f $temp_file

echo "  "
echo "  "
echo "Processing complete. Summary saved to $output_summary_file."
echo "  "
