form mat_amp.frm > withoutcolor.m && sed '1,4d' withoutcolor.m|sed 's/mat/L   mat/g' > colorform.m 

form colorfactor.frm > out.m

cat out.m

echo $'************************************************************ ' 
echo $' if out.m shows error check withoutcolor.m  ' :
echo $'************************************************************ ' 
#rm colorform.m withoutcolor.m
