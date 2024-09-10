#!/bin/bash

#script of extracting integrals from the amplitude for reduction 

rm integral.m  &&  rm integralsto.reduze 

for i in {1..10}
do
        sed -i '5s/.*/#$diaS = '$i';/g' input.h && sed -i '6s/.*/#$dia = '$i';/g' input.h

        tform -w10 -l mat_amp.frm
	bash extractintegral.sh

        echo "mat"$i" done!"    
done
