c---------------------------------------------------------------------
c       This file contains      
c     1.  Born / reduced_Born DrellYan
c     2.  1real DrellYan ( qq ,qg ) channels 
c---------------------------------------------------------------------
       subroutine uu2ee_r(p,msq)
       implicit double precision (a-h,o-z)
       dimension p(0:3,1:5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       parameter(Pi=3.141592653589793238D0)
       double precision msq(1:5),msq1,msq2
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
 
       msq(1)= ((-16*gs**2*qe**4*qu**2*(s14*(-2*s13*s25 + s23*(s15 +
     .        s25) - s15*s35) +
     .  s24*(-2*s15*s23 + s13*(s15 + s25) - s25*s35) -
     . (s13*s15 + s23*s25)*s45 + s12*(-2*s14*s23 - 2*s13*s24 + s14*s35 +
     .    s24*s35 + (s13 + s23)*s45)))/(9*s15*s25*s34**2))

c           msq2 = 1.0d0*(16*gs**2*qe**4*qu**2*
c     -    (s14**3*(-(s35*(s34 + 2*s35)) + s34*s45) +
c     -      s13**3*(s34*(s35 - s45) - 2*s45**2) +
c     -      s13**2*(s34*(s15*s34 + 3*s14*s35 + 2*s15*s35 - 2*s34*s35) +
c     -         s34*(s14 + 2*s34)*s45 + 4*s14*s35*s45 -
c     -         2*(s14 + s15 - 3*s34 - s35)*s45**2 + 2*s45**3) +
c     -      s15*s34**2*(2*s34**2 + s35**2 + 4*s35*s45 + s45**2 +
c     -         4*s34*(s35 + s45) - 2*s15*(s34 + s35 + s45)) +
c     -      s14**2*(6*s34*s35**2 + 2*s34**2*(s35 - s45) +
c     -         2*s35**2*(s35 + s45) +
c     -         s15*(s34**2 - 2*s35**2 + 2*s34*s45)) +
c     -      s14*s34*(-(s35**2*(2*s34 + s35)) +
c     -         (2*s34**2 + 4*s34*s35 + s35**2)*s45 +
c     -         (2*s34 + 3*s35)*s45**2 + s45**3 +
c     -         2*s15**2*(s34 + s35 + s45) -
c     -         2*s15*(2*s34**2 + 3*s34*s35 + 4*s34*s45 + 3*s35*s45 +
c     -            s45**2)) +
c     -      s13*(s14**2*(-2*s35*(s35 - 2*s45) + s34*(s35 + 3*s45)) +
c     -         s34*(2*s34**2*s35 + 2*s15**2*(s34 + s35 + s45) -
c     -            2*s15*(2*s34**2 + 4*s34*s35 + s35**2 +
c     -               3*(s34 + s35)*s45) +
c     -            2*s34*(s35**2 + 2*s35*s45 - s45**2) +
c     -            (s35 + s45)*(s35**2 + 2*s35*s45 - s45**2)) +
c     -         2*s14*(s15*(2*s34**2 + 2*s35*s45 + 3*s34*(s35 + s45)) -
c     -            2*(3*s34*s35*s45 + s34**2*(s35 + s45) +
c     -               s35*s45*(s35 + s45))))))/
c     -  (9.*(s13 + s14 - s34)**2*s34**2*(s34 + s35 + s45)**2)
           
c          msq2 =  (16*gs**2*qe**4*qu*
c     -    (s13**3*(s34*(s35 - s45) - 2*s45**2) + 
c     -      s14**3*(-2*s35**2 + s34*(-s35 + s45)) + 
c     -      s15*s34**2*(2*s34**2 + s35**2 + 4*s35*s45 + s45**2 + 
c     -         4*s34*(s35 + s45) - 2*s15*(s34 + s35 + s45)) + 
c     -      s14**2*(6*s34*s35**2 + 2*s34**2*(s35 - s45) + 
c     -         2*s35**2*(s35 + s45) + 
c     -         s15*(s34**2 - 2*s35**2 + 2*s34*s45)) + 
c     -      s13**2*(6*s34*s45**2 + 2*s34**2*(-s35 + s45) + 
c     -         2*s45**2*(s35 + s45) + 
c     -         s15*(s34**2 + 2*s34*s35 - 2*s45**2) + 
c     -         s14*(3*s34*s35 + s34*s45 + 4*s35*s45 - 2*s45**2)) + 
c     -      s14*s34*(-s35**3 + 2*s34**2*s45 + s35**2*s45 + 
c     -         3*s35*s45**2 + s45**3 + 2*s15**2*(s34 + s35 + s45) - 
c     -         2*s15*(2*s34**2 + 3*s34*s35 + 4*s34*s45 + 3*s35*s45 + 
c     -            s45**2) + s34*(-2*s35**2 + 4*s35*s45 + 2*s45**2)) + 
c     -      s13*(s14**2*(-2*s35*(s35 - 2*s45) + s34*(s35 + 3*s45)) + 
c     -         s34*(2*s34**2*s35 + s35**3 + 3*s35**2*s45 + s35*s45**2 - 
c     -            s45**3 + 2*s15**2*(s34 + s35 + s45) - 
c     -            2*s15*(2*s34**2 + 4*s34*s35 + s35**2 + 3*s34*s45 + 
c     -               3*s35*s45) + 2*s34*(s35**2 + 2*s35*s45 - s45**2))
c     -          + 2*s14*(s15*
c     -             (2*s34**2 + 2*s35*s45 + 3*s34*(s35 + s45)) - 
c     -            2*(3*s34*s35*s45 + s34**2*(s35 + s45) + 
c     -               s35*s45*(s35 + s45))))))/
c     -  (9.*(s13 + s14 - s34)**2*s34**2*(s34 + s35 + s45)**2)
c
         msq(2) = gs**2*qe**4*qu**2*((16*s14*s23)/(9.*s15*s34**2) + 
     -    (16*s13*s24)/(9.*s15*s34**2) + (32*s23*s24)/(9.*s12*s34**2) - 
     -    (16*s14*s23*s25)/(9.*s12*s15*s34**2) - 
     -    (16*s13*s24*s25)/(9.*s12*s15*s34**2) + 
     -    (16*s14*s35)/(9.*s12*s34**2) + (16*s24*s35)/(9.*s12*s34**2) - 
     -    (16*s24*s35)/(9.*s15*s34**2) + 
     -    (16*s14*s25*s35)/(9.*s12*s15*s34**2) + 
     -    (32*s24*s25*s35)/(9.*s12*s15*s34**2) + 
     -    (16*s13*s45)/(9.*s12*s34**2) + (16*s23*s45)/(9.*s12*s34**2) - 
     -    (16*s23*s45)/(9.*s15*s34**2) + 
     -    (16*s13*s25*s45)/(9.*s12*s15*s34**2) + 
     -    (32*s23*s25*s45)/(9.*s12*s15*s34**2) - 
     -    (32*s35*s45)/(9.*s15*s34**2))

         msq(3)= gs**2*qe**4*qu**2*((32*s13*s14)/(9.*s12*s34**2) +
     -    (16*s14*s23)/(9.*s25*s34**2) -
     -    (16*s14*s15*s23)/(9.*s12*s25*s34**2) +
     -    (16*s13*s24)/(9.*s25*s34**2) -
     -    (16*s13*s15*s24)/(9.*s12*s25*s34**2) +
     -    (16*s14*s35)/(9.*s12*s34**2) + (16*s24*s35)/(9.*s12*s34**2) -
     -    (16*s14*s35)/(9.*s25*s34**2) +
     -    (32*s14*s15*s35)/(9.*s12*s25*s34**2) +
     -    (16*s15*s24*s35)/(9.*s12*s25*s34**2) +
     -    (16*s13*s45)/(9.*s12*s34**2) + (16*s23*s45)/(9.*s12*s34**2) -
     -    (16*s13*s45)/(9.*s25*s34**2) +
     -    (32*s13*s15*s45)/(9.*s12*s25*s34**2) +
     -    (16*s15*s23*s45)/(9.*s12*s25*s34**2) -
     -    (32*s35*s45)/(9.*s25*s34**2))

            ! Total contribution from qg ang qq initiated process
       ! + msq2 !this msq2 is the contribution from qg initiaied process
c      msq = msq1

      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c     [u U -> e E]  Born 
c--------------------------------------------------------------------o
       function Born_uU2eE(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       common/usedalpha/AL,ge
c       ge=1d0/128d0
       e= DSQRT(ge*4.d0*PI)
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2
      s13 =  2.0d0*dot(p1,p3) ! t
      s23 =  2.0d0*dot(p1,p4) ! u
      s12 =  2.0d0*dot(p1,p2) ! s
      qu2 = 1d0!4d0/9d0

      t =  -2.0d0*dot(p1,p3) ! t
      u =  -2.0d0*dot(p1,p4) ! u
      s =  2.0d0*dot(p1,p2) ! s

c	s = 1000d0
c	t = -10000d0
c	u = -15000d0



c      Born_uU2eE= CF*(2*e**4*qu2*(-2*s13*s23 + s12*(s13 +
c     .            s23)))/(3d0*s12**2)
	Born_uU2eE =1* e**4*qu2*(t*t + u*u)/s/s/3d0
c	Born_uU2eE = e**4*qu2*(t/u + u/t)/3d0
c	write(*,*)'sig = ',Born_uU2eE
c	stop
       return
       end
c---------------------------------------------------------------------
