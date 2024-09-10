#procedure gammaprop

repeat,id G(?a,li1?,li2?,li3?,li4?,li1?,?b)=-2*G(?a,li4,li3,li2,?b)-(n-4)*G(?a,li2,li3,li4,?b); 
.sort
repeat,id G(?a,li1?,li2?,li3?,li1?,?b)=4*d_(li2,li3)*G(?a,?b)+(n-4)*G(?a,li2,li3,?b); 
.sort
repeat,id G(?a,li1?,li2?,li1?,?b)=(-n+2)*G(?a,li2,?b); 
.sort

#endprocedure
