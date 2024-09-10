#-
off statistics,finalstats,allwarnings;
nwrite statistics;

Symbol ch,mh,NA,n,AS,PI,v,s12,AL,flag;
Index mu,nu,a,b,c,d,al,bt;
Index
#do ii=1,20;
cix`ii'=NA,cix10`ii'=NA,lix`ii'=n,lix10`ii'=n
#enddo
      ;
Vector p1,p2,p3,p4,p5;
Index ca=NA,cb=NA,cal=NA,cbt=NA;
Vector p1,p2,q,p,nv;
CFunction epol,db;

Local   A = -i_*ch*db(cix1,p1)*db(cix3,p2)*(-d_(lix1,lix3)*p1.p2 + p1(lix3)*p2(lix1))*epol(lix1,p1)*epol(lix3,p2); 
Local  Ac = +i_*ch*db(cix101,p1)*db(cix103,p2)*(-d_(lix101,lix103)*p1.p2 + p1(lix103)*p2(lix101))*epol(lix101,p1)*epol(lix103,p2); 
*Local Ac= -ch*db(cal,p)*db(cbt,q)*(d_(al,bt)*mh^2/2 - p(bt)*q(al))*epol(al,p)*epol(bt,q); 
.sort

Local A2 = Ac*A*d_(cix1,cix3)*d_(cix101,cix103);
.sort
print,A2;
.sort
*id epol(mu?, p?) * epol(b?, p?) = -d_(mu, b) + (p(mu) * nv(b) + nv(b) * p(mu)) / p.nv;
id epol(lix1?,p1?)*epol(lix101?,p1?) = - d_(lix1,lix101)+flag*(p1(lix1)*nv(lix101)+nv(lix1)*p1(lix101))/p1.nv;
*id epol(lix1?,p1?)*epol(lix101?,p1?) = - d_(lix1,lix101) + (p1(lix1)*nv(lix101)+nv(lix101)*p1(lix1))/p1.nv;
*id epol(mu?,p?)*epol(b?,p?) = - d_(mu,b)+(p(mu)*nv(b)+nv(b)*p(mu))/p.nv;
.sort
print,A2;
.sort
id p1.p1 = 0;
id p2.p2 = 0;
id p1.p2 = s12/2;
id q.q = 0;
id nv.nv = 0;
id n=4;
*id ch = AS/3/PI/v;
.sort
B flag;
print +s,A2;
.sort
*id db(ca?,p?)*db(cal?,p?)=d_(ca,cal);
id db(cix1?,p1?)*db(cix101?,p1?)=d_(cix1,cix101);
.sort
*id NA=8;
id AS = AL/4/PI;
.sort
print +s,A2;
.sort

.end
