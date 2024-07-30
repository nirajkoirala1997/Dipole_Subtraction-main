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

id me = massofelectron;
id mu = massofup;
id gew=0;
.sort

id epol(lix1?,p1?,0)*epol(lix101?,p1?,0) = -d_(lix1,lix101)+(p1(lix1)*nv(lix101)+nv(lix101)*p1(lix1))/p1.nv;
.sort

*l pol1=-d_(lix6,lix106)+(p5(lix6)*nv(lix106)+nv(lix106)*p5(lix6))/p5.nv;

*multiply(pol1);

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
L coloconn =d_(cifx3,cifx103)*d_(cifx1,cifx101)*d_(cix6,cix106);

multiply coloconn;
repeat, id T(cifx1?,cifx2?,?a)*T(cifx2?,cifx3?,?b)=T(cifx1,cifx3,?a,?b);
id T(cifx1?,cifx1?,?a)=Tr(?a);
*********************************************************
************************
* MANDELSTAM VARIABLES *
************************
id fprop(p1-p4-p3)=1/(-s13-s14+s34);
id fprop(-p2+p3+p4)=1/(-s23-s24+s34);
id phprop(p3+p4)=1/(s34);
id fprop(p1-p5)=1/s15;
id fprop(-p2+p5)=1/s25;


#call mandelsterm2
************************
B fprop;
print +s mat;
#printtimes;
.end

