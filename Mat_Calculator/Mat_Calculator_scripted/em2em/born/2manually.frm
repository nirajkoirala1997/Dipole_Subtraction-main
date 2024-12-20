#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include ../../main_files/def.h
#include input.h
#include ../../main_files/feyn.h
#include mandelsterm.h
#include ../../main_files/grfunc.h
#include ../../main_files/SOn.prc
#include ../../main_files/SUn.prc
#include ../../main_files/color.h
#include ../../main_files/gamma5.h
#include ../../main_files/Camplitude.h
#include ../../main_files/amplitude.h
.sort
l mat  = amp*ampc;
.sort
id U(si1?,p1?,x2?)*UB(si2?,p1?,x2?)=G(si1,si2,p1)+x2*G(si1,si2);
id V(si1?,p1?,x2?)*VB(si2?,p1?,x2?)=G(si1,si2,p1)-x2*G(si1,si2);
.sort
repeat,id G(si1?,si2?,?a)*G(si2?,si3?,?b)=G(si1,si3,?a,?b);
.sort
#call mass
.sort
#do i=1,20;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo
.sort

*multiply replace_(p4,p1+p2-p3);
.sort
#do i=1,20;
tracen `i';
#enddo
.sort
id phprop(-p1 + p3) =1/t;
#call mandelsterm
.sort
B qe;
print +s mat;
#write <out.m> "      mat = %e",mat
.end


*************************
* MOMENTUM CONSERVATION *
*************************
multiply replace_(p4,p1+p2-p3);

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
.sort
.sort

id n = 4 ;
*id mH = 0;
*id mH^2 = 0;

#call mandelsterm
id p1.p1 = 0;
id p2.p2 = 0;
id p3.p3 = 0;
id p4.p4 = 0;
id gprop( - p1 - p2) = 1/s^2;
id fprop( - p1 + p3) = -1/t^2;
id fprop( - p2 + p3) = -1/u^2;
id fprop( + p2 - p4) = -1/t^2;
id fprop( + p1 - p4) = -1/u^2;
id phprop( -p1 + p3) = -1/t;
.sort
#include ../../main_files/colordef.h
#call SUn
.sort 
id a =1/2;
id nf=1;

id NA=NF^2-1;
id NF=3;
id 1/NF=1/3;
.sort
*id flag1  = 0;
*Format mathematica;
*B kg,p1.nv,p2.nv,p3.nv;
B qe;
print +s mat;
*#write <out.m> "%O"
#write <out.m> "      mat = %e",mat
#printtimes;
.end

