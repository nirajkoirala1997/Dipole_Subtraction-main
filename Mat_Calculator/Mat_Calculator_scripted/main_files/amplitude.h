
#do i=`$diaS',`$dia'
#include amp.qgraf #d`i'

l amp`i' = Rq ;
.sort
Drop Rq;
.sort
#enddo


Local amp=
#do i=`$diaS',`$dia'
+amp`i'
#enddo
;

#call feynrules

repeat,id G(si1?,si2?,?a)*G(si2?,si3?,?b)=G(si1,si3,?a,?b);
*repeat,id T(cif1?,cif2?,?a)*T(cif2?,cif3?,?b)=T(cif1,cif3,?a,?b);

