#-
off statistics,finalstats,allwarnings;
nwrite statistics;

#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/def.h
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/gamma5.h
#include /home/vaibhav/work/matrix_amp/NEWcode/main_files/gamma_prop.prc


L exp1 = G(six1,six1,p1,lix102,p1,lix104,p2,lix4,g5,p2,li8,k1,li8,p1,
      lix2);

repeat id G(?a,g5,li1?,?c,g5,?d,g5,?e,g5,?b) = -G(?a,li1,g5,?c,g5,?d,g5,?e,g5,?b);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c,g5,?d,g5,?b) = -G(?a,li1,g5,?c,g5,?d,g5,?b);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c,g5,?d) = -G(?a,li1,g5,?c,g5,?d);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c) = -G(?a,li1,g5,?c);

P ;
.end

