#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/def.h
#include input.h
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/feyn.h
#include mandelsterm.h
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/SOn.prc
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/SUn.prc
#include /home/vaibhav/work/matrix_amp/NEWcode2/main_files/kpositive.prc
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
id epolglu(lix1?,p1?,0)*epolglu(lix101?,p1?,0) = - d_(lix1,lix101)+(p1(lix1)*nv(lix101)+nv(lix101)*p1(lix1))/p1.nv;
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
*multiply replace_(p4,p1+p2-p3-p5);

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


#call mandelsterm
******************
* LOOP wala kaam 
******************

#include ./reduze/qQ2ph1L.out

id gprop(?a) = Prop(?a,0);
id fprop(?a) = Prop(?a,0);

#call DiaMatch`$diaS'
#call mandelsterm

***
*id Prop(-k1,0) = Prop(k1,0);
*id 1/Prop(-k1,0) = 1/Prop(k1,0);
**
*id Prop(-k1+p1,0) = Prop(k1 - p1,0);
*id 1/Prop(-k1+p1,0) = 1/Prop(k1 - p1,0);
**
*id Prop( - k1 + p1 + p2,0) = Prop(  k1 - p1 - p2,0);
*id 1/Prop( - k1 + p1 + p2,0) = 1/Prop(  k1 - p1 - p2,0);
**
**
*id Prop(k1,0) = propA;
*id 1/Prop(k1,0) = 1/propA;
**
*id Prop(k1 - p1,0) = propB;
*id 1/Prop(k1 - p1,0) = 1/propB;
**
*id Prop(  k1 - p1 - p2,0) = propC;
*id 1/Prop(  k1 - p1 - p2,0) = 1/propC;
**
*
id F1 = INT(F1);
id INT(?a)*propA^xx?=INT(?a,xx);
id INT(?a)*propB^xx?=INT(?a,xx);
id INT(?a)*propC^xx?=INT(?a,xx);
**
#include ./reduze/to_amp.frm


multiply replace_(d,n);

************************
B fprop;
print +s mat;
#printtimes;
.end

