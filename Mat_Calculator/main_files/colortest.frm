#-
off statistics,finalstats,allwarnings;
nwrite statistics;

*#include def.h
#include colordef.h
#include color.h
#include SUn.prc


*L exp2 = Tr(ci2,ci102)*f(ci2,cix102,cix104)*f(ci102,cix102,cix104);
*L exp2 = f(ci1,ci2,ci3);
L exp2 = Tr(ci1,ci2,ci3)*f(ci3,ci2,ci1)*i_;
p +s exp2;
.sort
#call SUn 

id NF=3;
id a =1/2;
id nf=1;
p +s exp2;
.sort

