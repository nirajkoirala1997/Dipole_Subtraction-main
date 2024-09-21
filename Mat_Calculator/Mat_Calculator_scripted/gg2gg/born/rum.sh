#*form mat_amp.frm > withoutcolor.m && sed '1,3d' withoutcolor.m|sed 's/mat/L   mat/g' > colorform.m 
form mat_amp.frm && sed '1,3d' withoutcolor.m|sed 's/mat/L   mat/g' > colorform.m 

form colorfactor.frm > out2.m 

#rm colorform.m
