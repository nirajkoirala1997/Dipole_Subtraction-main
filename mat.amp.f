c---------------------------------------------------------------------
       subroutine uu2ee_r(p,msq2)
       implicit double precision (a-h,o-z)
       dimension p(0:3,1:5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       parameter(Pi=3.141592653589793238D0)
       double precision msq2
       call p2dtop1d_5(p,p1,p2,p3,p4,p5)

       s12=2.d0*dot(p1,p2)
       s13=2.d0*dot(p1,p3)
       s14=2.d0*dot(p1,p4)
       s15=2.d0*dot(p1,p5)
       s23=2.d0*dot(p2,p3)
       s24=2.d0*dot(p2,p4)
       s25=2.d0*dot(p2,p5)
       s34=2.d0*dot(p3,p4)
       s35=2.d0*dot(p3,p5)
       s45=2.d0*dot(p4,p5)

       AL=0.118d0
       ge=0.007547169811320755d0
       qe=DSQRT(ge*4.d0*PI)
       gs=DSQRT(AL*4.d0*PI)
c       gs=AL
       qu = 2d0/3d0
 
      msq2= ((-16*gs**2*qe**4*qu**2*(s14*(-2*s13*s25 + s23*(s15 +
     .        s25) - s15*s35) +
     .  s24*(-2*s15*s23 + s13*(s15 + s25) - s25*s35) -
     . (s13*s15 + s23*s25)*s45 + s12*(-2*s14*s23 - 2*s13*s24 + s14*s35 +
     .    s24*s35 + (s13 + s23)*s45)))/(9*s15*s25*s34**2))
      return
      end
c---------------------------------------------------------------------
