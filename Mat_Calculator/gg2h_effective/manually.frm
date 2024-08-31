#-
off statistics,finalstats,allwarnings;
nwrite statistics;

Index mu,nu,a,b,c,d,al,bt,ca,cb,cal,cbt;
Vector p,q,nv;
CFunction epol,db;
Symbol ch,mh;

Local  A= -ch*db(ca,p)*db(cb,q)*(d_(mu,nu)*mh^2/2 - p(nu)*q(mu))*epol(mu,p)*epol(nu,q); 
Local Ac= -ch*db(cal,p)*db(cbt,q)*(d_(al,bt)*mh^2/2 - p(bt)*q(al))*epol(al,p)*epol(bt,q); 
.sort
Local A2 = Ac*A;
.sort
id epol(mu?,p?)*epol(b?,p?) = - d_(mu,b)+(p(mu)*nv(b)+nv(b)*p(mu))/p.nv;
.sort
print,A2;
.sort
id p.p = 0;
id q.q = 0;
id nv.nv = 0;
.sort
print +s,A2;
.sort
id db(ca?,p?)*db(cal?,p?)=d_(ca,cal);
.sort
print +s,A2;

.end
