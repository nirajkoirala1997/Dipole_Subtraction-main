#!/bin/bash

#script of extracting integrals from the amplitude for reduction 

rm integral.m  &&  rm integralsto.reduze 

for i in {1..1}
do
        sed -i '5s/.*/#$diaS = '$i';/g' input.h && sed -i '6s/.*/#$dia = '$i';/g' input.h

        tform -w10 -l mat_amp.frm > tmp.tmp
	grep -Po 'INT.*\)' mat_amp.log | sed 's/INT(//g' | sed 's/)//g' | sed 's/,/ /g' | sort | uniq >> integral.m

	grep  F.* integral.m | sort | uniq > integralsto.reduze

	rm tmp.tmp
        echo "mat"$i" done!"    
done

rm integral.m
