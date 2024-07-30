#!/bin/bash
#script form generating matrix amplitudes for all diagrams
#Here loop is running over number of diagrams for the amplitude only not for the its complex conjugate

rm amplitudes/mat*

for i in {1..10}
do 
	sed -i '9s/.*/#$diaS = '$i';/g' input.h && sed -i '10s/.*/#$dia = '$i';/g' input.h

	tform -w10 -l mat_amp.frm > amplitudes/mat"$i".m
	sed -i 's/mat/mat'$i'/g' amplitudes/mat"$i".m
        sed -i 's/=/=(/g' amplitudes/mat"$i".m
        sed -i 's/;/);/g' amplitudes/mat"$i".m
	
	echo "mat"$i" done!"	
echo "----------------------------------------------------------------- " 
done

sed -i '1,3d' amplitudes/mat*
echo "Baaki kaam khud karo !!,
      (ㆆ_ㆆ) ᕙ(^▿^-ᕙ) ( ◡́.◡̀)(^◡^ ) (◍•ᴗ•◍)"	
echo "----------------------------------------------------------------- " 
