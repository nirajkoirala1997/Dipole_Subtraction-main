#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/def.h
#include input.h
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/feyn.h
#include mandelsterm.h
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/SOn.prc
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/SUn.prc
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/color.h
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/gamma5.h
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/Camplitude.h
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/amplitude.h
.sort


l mat  = amp*ampc;


id U(si1?,p1?,x2?)*UB(si2?,p1?,x2?)=G(si1,si2,p1)+x2*G(si1,si2);
id V(si1?,p1?,x2?)*VB(si2?,p1?,x2?)=G(si1,si2,p1)-x2*G(si1,si2);
repeat,id G(si1?,si2?,?a)*G(si2?,si3?,?b)=G(si1,si3,?a,?b);
.sort

id gew=0;
id mu=0;
.sort

l pol1=-d_(lix3,lix103)+(p1(lix3)*nv(lix103)+nv(lix103)*p1(lix3))/p1.nv;

multiply(pol1);

#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo

#do i=1,20;
tracen `i';
#enddo

.sort


*************************
* MOMENTUM CONSERVATION *
*************************
*multiply replace_(p5,p1+p2-p3-p4);
*multiply replace_(p2,-p1+p3+p4+p5);

************
* SPIN AVG *
************
multiply (1/4);

*************
* COLOR AVG *
*************
multiply (1/9);

***************
* COLOR TRACE *
***************
.sort
L coloconn =d_(cix3,cix103)*d_(cifx1,cifx101)*d_(cifx6,cifx106);

multiply coloconn;
repeat, id T(cifx1?,cifx2?,?a)*T(cifx2?,cifx3?,?b)=T(cifx1,cifx3,?a,?b);
id T(cifx1?,cifx1?,?a)=Tr(?a);

*********************************************************
************************
* MANDELSTAM VARIABLES *
************************
id fprop(p1 + p2)=1/s12;
id phprop(p3 + p4)=1/s34;
id fprop( - p1 + p5)=-1/s15;
id fprop(p3 + p4 + p5)=1/(s34+s45+s35);
*id fprop(p2 - p3 - p4) = 1/(-s23-s24+s34);
#call mandelsterm2
************************
id n=4;
B fprop,phprop;
print +s mat;
#printtimes;
.end

