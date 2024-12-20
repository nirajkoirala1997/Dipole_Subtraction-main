#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include ../../main_files/def.h
.sort
Local amp =
      epolglu(lix1,p1,0)*epolglu(lix3,p2,0)*Vx(1,3,glugluHig,-1,-3,-2,p1,p2,
      -p3)*db(cix1,p1)*db(cix3,p2);
.sort
#do i = 1,3
#do j = 1,3
id Vx(x1?,x2?,glugluHig,-`i',-`j',x3?,p1?,p2?,p3?)  = -i_*ch*d_(cix`i',cix`j')*(-d_(lix`i',lix`j')*p1.p2 + p1(lix`j')*p2(lix`i'));
#enddo
#enddo
.sort

Local ampc =
      epolglu(lix101,p1,0)*epolglu(lix103,p2,0)*Vx(1,3,glugluHig,-1,-3,-2,p1,p2,
      -p3)*db(cix101,p1)*db(cix103,p2);
.sort
#do i = 1,3
#do j = 1,3
id Vx(x1?,x2?,glugluHig,-`i',-`j',x3?,p1?,p2?,p3?)  = +i_*ch*d_(cix10`j',cix10`i')*(-d_(lix10`j',lix10`i')*p2.p1 + p1(lix10`i')*p2(lix10`j'));
*id Vx(x1?,x2?,glugluHig,-`i',-`j',x3?,p1?,p2?,p3?)  = +i_*ch*d_(cix10`i',cix10`j')*(-d_(lix10`i',lix10`j')*p1.p2 + p1(lix10`j')*p2(lix10`i'));
#enddo
#enddo
.sort
Local mat = ampc*amp;
.sort
print, amp,ampc,mat;
.sort
id epolglu(lix1?,p1?,0)*epolglu(lix101?,p1?,0) = - d_(lix1,lix101)+(p1(lix1)*nv(lix101)+nv(lix101)*p1(lix1))/p1.nv;
.sort
print,mat;
.sort
id db(cix1?,p1?)*db(cix101?,p1?)=d_(cix1,cix101);

.sort
id nv.nv =0;
id p1.p1 =0;
id p2.p2 =0;
id p1.p2 =mH^2/2;
id n = 4;
.sort
print,+s mat;



.end















#call feynrules(amp)
#call feynrules(ampc)
.sort
print,amp,ampc;
.end
id in(glu(-1,p1)) = epolglu(lix1,p1,0)*db(cix1,p1);
id in(glu(-3,p2)) = epolglu(lix3,p2,0)*db(cix3,p2);

.end



id ou(Hig(-2,p1?)) = 1;
print,Rq;
.end


#include ../../main_files/feyn.h
#include mandelsterm.h
#include ../../main_files/grfunc.h
#include ../../main_files/SOn.prc
#include ../../main_files/SUn.prc
*#include ../../main_files/color.h
#include ../../main_files/gamma5.h
#include ../../main_files/Camplitude.h
#include ../../main_files/amplitude.h
.sort

l mat  = amp*ampc;
.sort
print,amp,ampc;
.end


*id U(si1?,p1?,x2?)*UB(si2?,p1?,x2?)=G(si1,si2,p1)+x2*G(si1,si2);
*id V(si1?,p1?,x2?)*VB(si2?,p1?,x2?)=G(si1,si2,p1)-x2*G(si1,si2);
*repeat,id G(si1?,si2?,?a)*G(si2?,si3?,?b)=G(si1,si3,?a,?b);
*.sort
*id mgr = massofgraviton;
*#call mass

*#call grfunc
.sort
*id epolph(lix1?,p3?,0)*epolph(lix101?,p3?,0) = - d_(lix1,lix101);
id epolglu(lix1?,p1?,0)*epolglu(lix101?,p1?,0) = - d_(lix1,lix101)+(p1(lix1)*nv(lix101)+nv(lix101)*p1(lix1))/p1.nv;
.sort
*print,mat;

*#do i=1,10;
*id,once G(six1?,six1?,?a)=g_(`i',?a) ;
*#enddo
.sort

*#do i=1,20;
*tracen `i';
*#enddo

.sort


*************************
* MOMENTUM CONSERVATION *
*************************
multiply replace_(p3,p1+p2);

***************
* COLOR TRACE *
***************
.sort
print,mat;
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
print,mat;
.end

id n = 4 ;

#call mandelsterm2
id p1.p1 = 0;
id p2.p2 = 0;

.sort
print, mat;
.sort
*id flag1  = 1;
.sort
*.sort
*Format mathematica;
*B kg,p1.nv,p2.nv,p3.nv;
**print +s mat;
*#write <out.m> "%O"
*#write <out.m> "      mat = %e",mat
#printtimes;
.end

