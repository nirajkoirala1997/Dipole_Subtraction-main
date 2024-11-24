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
*id mgr = massofgraviton;
#call mass

*#call grfunc
.sort
*id epolph(lix1?,p3?,0)*epolph(lix101?,p3?,0) = - d_(lix1,lix101);
id epolglu(lix1?,p1,0)*epolglu(lix101?,p1,0) = - d_(lix1,lix101)+flag1*(p1(lix1)*p2(lix101)+p1(lix101)*p2(lix1))/p1.p2;
id epolglu(lix1?,p2,0)*epolglu(lix101?,p2,0) = - d_(lix1,lix101)+flag2*(p2(lix1)*p1(lix101)+p2(lix101)*p1(lix1))/p2.p1;
id epolglu(lix1?,p4,0)*epolglu(lix101?,p4,0) = - d_(lix1,lix101)+flag3*(p4(lix1)*nv(lix101)+p4(lix101)*nv(lix1))/p4.nv;
.sort
*print +s mat;
*.end
id 1/p1.p2 = 1/p1p2;
id 1/p2.p1 = 1/p2p1;
id 1/p4.nv = 1/p4nv;
.sort
*#do i=1,20;
*id,once G(six1?,six1?,?a)=g_(`i',?a) ;
*#enddo
*.sort

*#do i=1,20;
*tracen `i';
*#enddo
*.sort


*************************
* MOMENTUM CONSERVATION *
*************************
multiply replace_(p3,p1+p2-p4);

***************
* COLOR TRACE *
***************
.sort
id df(cifx1?,p1?)*df(cifx101?,p1?)=d_(cifx1,cifx101);
id db(cix1?,p1?)*db(cix101?,p1?)=d_(cix1,cix101);
repeat, id T(cifx1?,cifx2?,?a)*T(cifx2?,cifx3?,?b)=T(cifx1,cifx3,?a,?b);
id T(cifx1?,cifx1?,?a)=Tr(?a);
.sort
*********************************************************
************************
* MANDELSTAM VARIABLES *
************************
id nv.nv=0;
.sort

id n = 4 ;
*id mH = 0;
*id mH^2 = 0


#call mandelsterm
id p1.p1 = 0;
id p2.p2 = 0;
id p3.p3 = 0;
id p4.p4 = 0;
*id gprop( - p1 - p2) = 1/s^2;
*id gprop( - p1 + p3) = -1/t^2;
*id gprop( + p2 - p3) = -1/u^2;
*id gprop( + p2 - p4) = -1/t^2;
*id gprop( + p1 - p4) = -1/u^2;
id gprop( - p1 - p2) = propm1m2;
id gprop( - p1 + p3) = propm13;
id gprop( + p2 - p3) = prop2m3;
id gprop( + p2 - p4) = prop2m4;
id gprop( + p1 - p4) = prop1m4;
id gprop( - p1 + p4) = propm14;
.sort

#include ../../main_files/colordef.h
#call SUn
.sort 
id a =1/2;
id nf=1;

id NA=NF^2-1;
id NF=3;
id 1/NF=1/3;
id flag1  = 1;
id flag2  = 1;
id flag3  = 1;
*id u = mH^2 - s - t;
*id 1/u = 1/(mH^2 - s - t);
*multiply (p4.nv);
*multiply (p2.nv);
*multiply (p1.nv);
.sort
*Format mathematica;
*B f, kg,p1.nv,p2.nv,p3.nv;
print +s mat;
*#write <out.m> "%O"
*#write <out.m> "      mat = %e",mat
#printtimes;
.end

