#!/bin/bash

one=1;
two=2;

fileinput=qgrafin.dat
echo "Tell the order"
read input

if [ "$input" == "1" ]
then
sed -e "s/loops= \([0-9]\);/loops= $input;/" -e "s/vsum\[gs, \([0-9]\), \([0-9]\)\];/vsum\[gs, 4, 4\];/" -e "s/\([a-zA-Z]\)\([0-9]\).qgraf/\1$input.qgraf/"  $fileinput > qgraf.dat 
fi

if [ "$input" == "2" ]
then
sed -e "s/loops= \([0-9]\);/loops= $input;/" -e "s/vsum\[gs, \([0-9]\), \([0-9]\)\];/vsum\[gs, 6, 6\];/" -e "s/\([a-zA-Z]\)\([0-9]\).qgraf/\1$input.qgraf/" $fileinput > qgraf.dat
fi
