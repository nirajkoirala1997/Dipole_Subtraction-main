#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/def.h
#include input.h
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/feyn.h
#include mandelsterm.h
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/SOn.prc
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/SUn.prc
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/color.h
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/gamma5.h
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/Camplitude.h
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/amplitude.h
.sort


l mat  = amp*ampc;
.sort

id U(si1?,p1?,x2?)*UB(si2?,p1?,x2?)=G(si1,si2,p1)+x2*G(si1,si2);
id V(si1?,p1?,x2?)*VB(si2?,p1?,x2?)=G(si1,si2,p1)-x2*G(si1,si2);
repeat,id G(si1?,si2?,?a)*G(si2?,si3?,?b)=G(si1,si3,?a,?b);
.sort

#call mass

id epolph(lix1?,p3?,0)*epolph(lix101?,p3?,0) = - d_(lix1,lix101);
id epolglu(lix1?,p1?,0)*epolglu(lix101?,p1?,0) = - d_(lix1,lix101)+(p1(lix1)*nv(lix101)+p1(lix1)*nv(lix101))/p1.nv;
.sort


#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo

#do i=1,20;
tracen `i';
#enddo

.sort

***************
*SPIN AVG* 
***************
multiply(1/4);
***************
*COLOR AVG* 
***************
multiply(1/9);
*************************
* MOMENTUM CONSERVATION *
*************************
multiply replace_(p4,p1+p2-p3-p5);

***************
* COLOR TRACE *
***************
.sort
id df(cifx1?,p1?)*df(cifx101?,p1?)=d_(cifx1,cifx101);
id db(cix1?,p1?)*db(cix101?,p1?)=d_(cix1,cix101);
repeat, id T(cifx1?,cifx2?,?a)*T(cifx2?,cifx3?,?b)=T(cifx1,cifx3,?a,?b);
id T(cifx1?,cifx1?,?a)=Tr(?a);
*********************************************************
************************
* MANDELSTAM VARIABLES *
************************
id nv.nv=0;


*id n = 4 ;

#call mandelsterm
************************
B fprop;
print +s mat;
#printtimes;
.end

