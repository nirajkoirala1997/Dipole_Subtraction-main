echo "Enter the name of the logfile:"
read filename
date > "$filename".log
echo $'\n<=  INPUT FILE  =>  \n' >> "$filename".log
cat input.input >> "$filename".log 
echo $'\n  <=  VEGAS INPUT  => \n' >> "$filename".log
cat vegas.input >> "$filename".log
echo $'\n  <=  DATA  => \n' >> "$filename".log
make && ./run.x | tee -a "$filename".log
date >> "$filename".log
mv "$filename".log logs/
