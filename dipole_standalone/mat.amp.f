c---------------------------------------------------------------------
       subroutine uu2ee_r(p,msq)
       implicit double precision (a-h,o-z)
       dimension p(0:3,1:5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       parameter(Pi=3.141592653589793238D0)
       double precision msq(1:2),msq1,msq2
       common/usedalpha/AL,ge
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

       qe=DSQRT(ge*4.d0*PI)
       gs=DSQRT(AL*4.d0*PI)
c       gs=AL
       qu =1d0! 2d0/3d0
 
       msq1= ((-16*gs**2*qe**4*qu**2*(s14*(-2*s13*s25 + s23*(s15 +
     .        s25) - s15*s35) +
     .  s24*(-2*s15*s23 + s13*(s15 + s25) - s25*s35) -
     . (s13*s15 + s23*s25)*s45 + s12*(-2*s14*s23 - 2*s13*s24 + s14*s35 +
     .    s24*s35 + (s13 + s23)*s45)))/(9*s15*s25*s34**2))

           msq2 = (16*gs**2*qe**4*qu**2*
     -    (s14**3*(-(s35*(s34 + 2*s35)) + s34*s45) +
     -      s13**3*(s34*(s35 - s45) - 2*s45**2) +
     -      s13**2*(s34*(s15*s34 + 3*s14*s35 + 2*s15*s35 - 2*s34*s35) +
     -         s34*(s14 + 2*s34)*s45 + 4*s14*s35*s45 -
     -         2*(s14 + s15 - 3*s34 - s35)*s45**2 + 2*s45**3) +
     -      s15*s34**2*(2*s34**2 + s35**2 + 4*s35*s45 + s45**2 +
     -         4*s34*(s35 + s45) - 2*s15*(s34 + s35 + s45)) +
     -      s14**2*(6*s34*s35**2 + 2*s34**2*(s35 - s45) +
     -         2*s35**2*(s35 + s45) +
     -         s15*(s34**2 - 2*s35**2 + 2*s34*s45)) +
     -      s14*s34*(-(s35**2*(2*s34 + s35)) +
     -         (2*s34**2 + 4*s34*s35 + s35**2)*s45 +
     -         (2*s34 + 3*s35)*s45**2 + s45**3 +
     -         2*s15**2*(s34 + s35 + s45) -
     -         2*s15*(2*s34**2 + 3*s34*s35 + 4*s34*s45 + 3*s35*s45 +
     -            s45**2)) +
     -      s13*(s14**2*(-2*s35*(s35 - 2*s45) + s34*(s35 + 3*s45)) +
     -         s34*(2*s34**2*s35 + 2*s15**2*(s34 + s35 + s45) -
     -            2*s15*(2*s34**2 + 4*s34*s35 + s35**2 +
     -               3*(s34 + s35)*s45) +
     -            2*s34*(s35**2 + 2*s35*s45 - s45**2) +
     -            (s35 + s45)*(s35**2 + 2*s35*s45 - s45**2)) +
     -         2*s14*(s15*(2*s34**2 + 2*s35*s45 + 3*s34*(s35 + s45)) -
     -            2*(3*s34*s35*s45 + s34**2*(s35 + s45) +
     -               s35*s45*(s35 + s45))))))/
     -  (9.*(s13 + s14 - s34)**2*s34**2*(s34 + s35 + s45)**2)
            ! Total contribution from qg ang qq initiated process
       msq(1) = msq1   ! + msq2 !this msq2 is the contribution from qg initiaied process
       msq(2) = msq2
c      msq = msq1

      return
      end
c---------------------------------------------------------------------
