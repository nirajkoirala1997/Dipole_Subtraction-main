#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include ../../../main_files/def.h
#include input.h
*#include ../../../main_files/feyn.h
#include ../../../main_files/fyncolorno.h
#include mandelsterm.h
#include ../../../main_files/grfunc.h
#include ../../../main_files/SOn.prc
#include ../../../main_files/SUn.prc
#include ../../../main_files/color.h
#include ../../../main_files/gamma5.h
#include ../../../main_files/Camplitude.h
#include ../../../main_files/amplitude.h
.sort
l mat  = amp*ampc;
.sort

#call mass
.sort

*id epolph(lix1?,p3?,0)*epolph(lix101?,p3?,0) = - d_(lix1,lix101);
id epolglu(lix1?,p1?,0)*epolglu(lix101?,p1?,0) = - d_(lix1,lix101)+flag1*(p1(lix1)*nv(lix101)+p1(lix101)*nv(lix1))/p1.nv;

multiply replace_(p4,p1+p2-p3);

id 1/p1.nv = 1/p1nv;
id 1/p2.nv = 1/p2nv;
id 1/p3.nv = 1/p3nv;
id 1/p4.nv = (p1.nv+p2.nv-p3.nv)/p4nv^2;
.sort

*************************
* MOMENTUM CONSERVATION *
*************************

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
.sort
id gprop( - p1 - p2) = 1/s;
id gprop( - p1 - p2) = 1/s;
id gprop( - p1 - p2) = 1/s;
id gprop( - p1 + p3) = 1/t;
id gprop(   p2 - p3) = 1/u;
.sort
*id flag1  = 0;
*Format mathematica;
*B kg,p1.nv,p2.nv,p3.nv;
*print +s mat;
*#write <out.m> "%O"
#write <out.m> "      mat = %e",mat
#printtimes;
.end

B gs,f,gprop;
print mat;
.end
