#procedure grfunc
id Bgr(li1?,li2?,li3?,li4?,k1?)=
                        (d_(li1,li3)-flag*k1(li1)*k1(li3)/s)*
                        (d_(li2,li4)-flag*k1(li2)*k1(li4)/s)
                       +(d_(li1,li4)-flag*k1(li1)*k1(li4)/s)*
                        (d_(li2,li3)-flag*k1(li2)*k1(li3)/s)
                   -2/3*(d_(li1,li2) -flag*k1(li1)*k1(li2)/s)*
                        (d_(li3,li4)-flag*k1(li3)*k1(li4)/s);
.sort
id Cgr(li1?,li2?,li3?,li4?)=
                         d_(li1,li3)*d_(li2,li4)
                        +d_(li1,li4)*d_(li2,li3)
                        -d_(li1,li2)*d_(li3,li4);

.sort
id Dgr(li1?,li2?,li3?,li4?,k1?,k2?)=
                  d_(li1,li2)*k1(li4)*k2(li3)
               - (d_(li1,li4)*k1(li2)*k2(li3)
                  +d_(li1,li3)*k1(li4)*k2(li2)
                  -d_(li3,li4)*k1(li1)*k2(li2)
                  +d_(li2,li4)*k1(li1)*k2(li3)
                  +d_(li2,li3)*k1(li4)*k2(li1)
                  -d_(li3,li4)*k1(li2)*k2(li1));

.sort
id Egr(li1?,li2?,li3?,li4?,k1?,k2?)=
                d_(li1,li2)*(k1(li3)*k1(li4)
               +k2(li3)*k2(li4)+k1(li3)*k2(li4))
               -(d_(li2,li4)*k1(li1)*k1(li3)+d_(li2,li3)*k2(li1)*k2(li4)
               +d_(li1,li4)*k1(li2)*k1(li3)+d_(li1,li3)*k2(li2)*k2(li4));
.sort
#endprocedure grfunc
