#procedure Coloring


id T(cif1?,cif2?,ci1?)*T(cif2?,cif3?,ci2?)*T(cif3?,cif1?,ci3?) = 1/4*(d(ci1,ci2,ci3)+i_*f(ci1,ci2,ci3));

id T(cif1?,cif2?,ci1?)*T(cif3?,cif4?,ci1?) = (1/2)*(d_(cif1,cif4)*d_(cif2,cif3)-(1/NF)*d_(cif1,cif2)*d_(cif3,cif4));

id T(cif1?,cif2?,ci1?)*T(cif2?,cif1?,ci2?) = (1/2)*d_(ci1,ci2);

id T(cif1?,cif1?,ci1?)=0;

id f(ci1?,ci2?,ci3?)*f(ci4?,ci2?,ci3?)=NF*d_(ci1,ci4);

id f(ci1?,ci4?,ci5?)*f(ci2?,ci5?,ci6?)*f(ci3?,ci6?,ci4?) = (NF/2)*f(ci1,ci2,ci3);

id f(ci1?,ci3?,ci4?)*d(ci2?,ci3?,ci4?)=0;

#endprocedure



#procedure coloring

repeat,id T(cif1?,cif2?,?a)*T(cif2?,cif3?,?b)=T(cif1,cif3,?a,?b);

id T(cif1?,cif1?,?a)=Tr(?a);

id f(ci1?,ci2?,ci3?)=-2*i_*(Tr(ci1,ci2,ci3)-Tr(ci1,ci3,ci2));

id d(ci1?,ci2?,ci3?)=2*(Tr(ci1,ci2,ci3)+Tr(ci1,ci3,ci2));


*id T(cifx1?,cifx2?,ci1?)*T(cifx2?,cifx3?,ci2?)*T(cifx3?,cifx4?,ci1?)*T(cifx4?,cifx1?,ci3?)=-(1/4/NF)*d_(ci2,ci3);

*id T(cifx1?,cifx2?,ci1?)*T(cifx2?,cifx3?,ci2?)*T(cifx3?,cifx1?,ci3?)=1/4*(d(ci1,ci2,ci3)+i_*f(ci1,ci2,ci3));

*id T(cifx1?,cifx2?,ci1?)*T(cifx3?,cifx4?,ci1?)=(1/2)*d_(cifx1,cifx4)*d_(cifx2,cifx3)-(1/NF)*d_(cifx1,cifx2)*d_(cifx3,cifx4);

*id T(cifx1?,cifx2?,ci1?)*T(cifx2?,cifx1?,ci2?)=(1/2)*d_(ci1,ci2);

*id T(cifx1?,cifx2?,ci1?)*T(cifx3?,cifx4?,ci1?)=(1/2)*(d_(cifx1,cifx4)*d_(cifx2,cifx3)-(1/NF)*(d_(cifx1,cifx2)*d_(cifx3,cifx4)));

*id Tr(ci1?,ci2?,ci2?,ci3?)=(NF/2-1/NF)*Tr(ci1,ci3);

*id Tr(ci1?,ci2?,ci1?,ci3?)=-(1/4/NF)*d_(ci2,ci3);

*id Tr(ci1?,ci2?,ci3?)=1/4*(d(ci1,ci2,ci3)+i_*f(ci1,ci2,ci3));

*id Tr(ci1?,ci2?)=(1/2)*d_(ci1,ci2);


****************************************************

id f(ci1?,ci2?,ci3?)*f(ci4?,ci2?,ci3?)=NF*d_(ci1,ci4);

id f(ci1?,ci2?,ci3?)*d(ci4?,ci2?,ci3?)=0;


#endprocedure
