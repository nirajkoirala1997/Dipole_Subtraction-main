form mat_amp.frm > withoutcolor.m && sed '1,4d' withoutcolor.m|sed 's/mat/L   mat/g' > colorform.m 

form colorfactor.frm > out.m

#rm colorform.m
