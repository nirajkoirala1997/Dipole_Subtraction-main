***********************
* SPACETIME DIMENSION *
***********************
Dimension n;

Symbol N,NA,NF,im,I,d;

****************
* LOOP MOMENTA *
****************
Vector k1,k2;

********************
* EXTERNAL MOMENTA *
********************
Vector p1,p2,p3,p4,p5,p6,p7,p8,p45,p35,l1,l2,q;

***************
* GAUGE CHECK *
***************
Vector nn,nv,nm;

*****************
* LORENTZ INDEX *
*****************
Index
#do ii=1,200;
 li`ii'=n,lix`ii'=n,
#enddo
      ;

**************
* SPIN INDEX *
**************
Index
#do ii=1,200;
 si`ii'=n,six`ii'=n,
#enddo
      ;

****************
* COLOR INDEX  *
****************
*********
* GLUON *
*********
Index 
#do ii=1,200;
ci`ii'=NA,cix`ii'=NA,
#enddo
      ;

*********
* QUARK *
*********
Index 
#do ii=1,200;
 cif`ii'=NF,cifx`ii'=NF,
#enddo
      ;

********************
* FOR COLOR TRACES *
********************
CFunction Tr(cyclic);
CFunction T, Tp, f(antisymmetric);

Index aa;
Index 
#do ii=1,200;
a`ii'=n,bb`ii'=n,cc`ii'=n,dd`ii'=n,ee`ii'=n,ff`ii'=n,gg`ii'=n,hh`ii'=n,ii`ii'=n,jj`ii'=n,kk`ii'=n,ll`ii'=n,mm`ii'=n,
#enddo
      ;
Index g5;

CFunction in,ou,upq,UPQ,elt,ELT,glu,ph,gr;
CFunction U,UB,V,VB,epolph,epolglu,epolmz,epolgr;
CFunction AA,Vx,zz,QQ,EE,GG,ggg,ggH,gggH;
Cfunction phprop,zprop,fprop,gprop,Prop,Iprop,grprop;
CFunction G,eps,G1,G2,G3,T;
CFunction db,df;
CFunction SP,J,INT,Den;

*
*

Symbol xx,x;
Symbol
#do ii=1,200;
x`ii',xx`ii',
#enddo
      ;
Symbol ELTeltzbos,UPQupqzbos,ELTeltph,UPQupqph,UPQupqglu,glugluglu,gluglugr,phphgr;
Symbol m,mz,me,mu,mH,mgr;
Symbol s,t,u;
Symbol qe,qu,gew,sw,cw,gs,ch,kg;
Symbol cf,ca;
Symbol cve,cae,cvu,cau;
*Symbol s12,s13,s14,s15;
*Symbol s21,s23,s24,s25;
*Symbol s31,s32,s34,s35;
*Symbol s41,s42,s43,s45;
*Symbol s51,s52,s53,s54;
*Symbol s11,s22,s33,s44,s55;
Symbol t12,t13,t14,t15;
Symbol t21,t23,t24,t25;
Symbol t31,t32,t34,t35;
Symbol t41,t42,t43,t45;
Symbol u12,u13,u14,u15;
Symbol u21,u23,u24,u25;
Symbol u31,u32,u34,u35;
Symbol u41,u42,u43,u45;
Symbol t11,t22,t33,t44,t55;
AutoDeclare Symbol s,t,u;
Symbol propA,propB,propC,propD;
Symbol k11,k12;
AutoDeclare Symbol F1,F1x12,F1x123,F1x124,F2,F3,F1x132,F1;
