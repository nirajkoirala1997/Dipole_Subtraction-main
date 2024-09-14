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
repeat,id G(si1?,si2?,?a)*G(si2?,si3?,?b)=G(si1,si3,?a,?b);
.sort
id mgr = massofgraviton;
*#call mass
.sort
#call grfunc
.sort
*print,mat;
*.end

*id epolglu(lix1?,p1?,0)*epolglu(lix101?,p1?,0) = - d_(lix1,lix101)+flag1*(p1(lix1)*nv(lix101)+p1(lix101)*nv(lix1))/p1.nv;
id epolph(lix1?,p3?,0)*epolph(lix101?,p3?,0) = - d_(lix1,lix101);
id epolglu(lix1?,p1,0)*epolglu(lix101?,p1,0) = - d_(lix1,lix101)+flag1*(p1(lix1)*n1(lix101)+p1(lix101)*n1(lix1))/p1.n1;
id epolglu(lix1?,p2,0)*epolglu(lix101?,p2,0) = - d_(lix1,lix101)+flag2*(p2(lix1)*n2(lix101)+p2(lix101)*n2(lix1))/p2.n2;
.sort
multiply 2*(p1.n1)/s; 
multiply 2*(p2.n2)/s; 
.sort
id n1 = p2;
id n2 = p1;
id flag = 1;
id flag1 = 1;
id flag2 = 1;
.sort
#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo
.sort

#do i=1,20;
tracen `i';
#enddo

.sort

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
*********************************************************
************************
* MANDELSTAM VARIABLES *
************************
.sort
id nv.nv=0;
id n1.n1=0;
id n2.n2=0;
id n1.n2=0;
id n2.n1=0;

id n = 4 ;
id NA = 8 ;

.sort
#call mandelsterm
.sort
*id phprop(- p1 -p2) = 1/s;
*id grprop( - p1 - p2) = 1/(s-mH^2);
id phprop(- p1 -p2) = 1;
id grprop(- p1 -p2) = 1;
.sort
************************
*id flag  = 1;
*print, mat;
*.end
*id flag1  = 1;
.sort
*id u = -s-t;
.sort
Format mathematica;
*B kg,p1.nv,p2.nv,p3.nv;
*print +s mat;
#write <out.m> "%O"
#write <out.m> "      mat = %e",mat
#printtimes;
.end

