c---------------------------------------------------------------------
c       Born and 1real Diphoton mat are here
c---------------------------------------------------------------------
       subroutine uu2eE_r(p,msq)
       implicit double precision (a-h,o-z)
       dimension p(0:3,1:5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       parameter(Pi=3.141592653589793238D0)
       double precision N,N2m1,msq(1:5),msq1,msq2
       common/usedalpha/AL,ge
       call p2dtop1d_5(p,p1,p2,p3,p4,p5)

       s12=2.d0*dot(p1,p2)
       s13=-2.d0*dot(p1,p3)
       s14=-2.d0*dot(p1,p4)
       s15=-2.d0*dot(p1,p5)
       s23=-2.d0*dot(p2,p3)
       s24=-2.d0*dot(p2,p4)
       s25=-2.d0*dot(p2,p5)
       s34=2.d0*dot(p3,p4)
       s35=2.d0*dot(p3,p5)
       s45=2.d0*dot(p4,p5)



       qe=DSQRT(ge*4.d0*PI)
       gs=DSQRT(AL*4.d0*PI)
       e = qe
       as = AL/4d0/Pi 
       N = 3d0
       Cf=(N*N-1.0D0)/2.D0/N
       Tf=1.0D0/2.0D0
       N2m1=N*N-1.0D0
c       gs=AL
       qu =1d0! 2d0/3d0
 
c       call SMQQB(S12,S13,S14,S15,S23,S24,S25,
c     &                  S34,S35,S45,SSMqqb)
c       msq1 = SSMqqb

        msq1 =  Pi**2*e**4*as*Cf/N * (  - 64.D0*S13**(-1)*S14**(-1)*
     &    S23**(-1)*S24**(-1)*S34**3 - 64.D0*S13**(-1)*S14**(-1)*
     &    S23**(-1)*S34**2 - 64.D0*S13**(-1)*S14**(-1)*S24**(-1)*S34**2
     &     + 32.D0*S13**(-1)*S14**(-1)*S15*S23**(-1)*S25**(-1)*S34**2
     &     - 64.D0*S13**(-1)*S14**(-1)*S15*S23**(-1)*S34 + 32.D0*
     &    S13**(-1)*S14**(-1)*S15*S24**(-1)*S25**(-1)*S34**2 - 64.D0*
     &    S13**(-1)*S14**(-1)*S15*S24**(-1)*S34 + 16.D0*S13**(-1)*
     &    S14**(-1)*S15*S25**(-1)*S34 - 16.D0*S13**(-1)*S14**(-1)*S15
     &     - 32.D0*S13**(-1)*S14**(-1)*S15**2*S23**(-1) - 32.D0*
     &    S13**(-1)*S14**(-1)*S15**2*S24**(-1) - 16.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S24**(-1)*S34**3 - 16.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S24**(-1)*S25*S34**2 + 64.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S25**(-1)*S34**3 - 144.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S34**2 + 80.D0*S13**(-1)*S15**(-1)*
     &    S23**(-1)*S25*S34 - 32.D0*S13**(-1)*S15**(-1)*S23**(-1)*
     &    S25**2 )
      msq1 = msq1 + Pi**2*e**4*as*Cf/N * (  - 16.D0*
     &    S13**(-1)*S15**(-1)*S24**(-1)*S34**2 - 16.D0*S13**(-1)*
     &    S15**(-1)*S24**(-1)*S25*S34 + 32.D0*S13**(-1)*S15**(-1)*
     &    S25**(-1)*S34**2 - 32.D0*S13**(-1)*S15**(-1)*S34 - 32.D0*
     &    S13**(-1)*S15**(-1)*S24*S25**(-1)*S34 + 32.D0*S13**(-1)*
     &    S15**(-1)*S24 - 48.D0*S13**(-1)*S23**(-1)*S24**(-1)*S34**2 -
     &    16.D0*S13**(-1)*S23**(-1)*S24**(-1)*S25*S34 - 32.D0*S13**(-1)
     &    *S23**(-1)*S24**(-1)*S25**2 - 128.D0*S13**(-1)*S23**(-1)*
     &    S25**(-1)*S34**2 + 64.D0*S13**(-1)*S23**(-1)*S34 - 64.D0*
     &    S13**(-1)*S23**(-1)*S25 - 64.D0*S13**(-1)*S24**(-1)*S34 - 32.D
     &    0*S13**(-1)*S24**(-1)*S25 - 64.D0*S13**(-1)*S25**(-1)*S34 +
     &    32.D0*S13**(-1) + 32.D0*S13**(-1)*S24*S25**(-1) + 128.D0*
     &    S13**(-1)*S15*S23**(-1)*S25**(-1)*S34 - 64.D0*S13**(-1)*S15*
     &    S23**(-1) + 32.D0*S13**(-1)*S15*S24**(-1)*S25**(-1)*S34 - 32.D
     &    0*S13**(-1)*S15*S24**(-1) + 48.D0*S13**(-1)*S15*S25**(-1) -
     &    32.D0*S13**(-1)*S15**2*S23**(-1)*S25**(-1) )
      msq1 = msq1 + Pi**2*e**4*as*Cf/N * (  - 16.D0*
     &    S14**(-1)*S15**(-1)*S23**(-1)*S24**(-1)*S34**3 - 16.D0*
     &    S14**(-1)*S15**(-1)*S23**(-1)*S24**(-1)*S25*S34**2 - 16.D0*
     &    S14**(-1)*S15**(-1)*S23**(-1)*S34**2 - 16.D0*S14**(-1)*
     &    S15**(-1)*S23**(-1)*S25*S34 + 64.D0*S14**(-1)*S15**(-1)*
     &    S24**(-1)*S25**(-1)*S34**3 - 144.D0*S14**(-1)*S15**(-1)*
     &    S24**(-1)*S34**2 + 80.D0*S14**(-1)*S15**(-1)*S24**(-1)*S25*
     &    S34 - 32.D0*S14**(-1)*S15**(-1)*S24**(-1)*S25**2 + 32.D0*
     &    S14**(-1)*S15**(-1)*S25**(-1)*S34**2 - 32.D0*S14**(-1)*
     &    S15**(-1)*S34 - 32.D0*S14**(-1)*S15**(-1)*S23*S25**(-1)*S34
     &     + 32.D0*S14**(-1)*S15**(-1)*S23 - 48.D0*S14**(-1)*S23**(-1)*
     &    S24**(-1)*S34**2 - 16.D0*S14**(-1)*S23**(-1)*S24**(-1)*S25*
     &    S34 - 32.D0*S14**(-1)*S23**(-1)*S24**(-1)*S25**2 - 64.D0*
     &    S14**(-1)*S23**(-1)*S34 - 32.D0*S14**(-1)*S23**(-1)*S25 - 128.
     &    D0*S14**(-1)*S24**(-1)*S25**(-1)*S34**2 + 64.D0*S14**(-1)*
     &    S24**(-1)*S34 )
      msq1 = msq1 + Pi**2*e**4*as*Cf/N * (  - 64.D0*
     &    S14**(-1)*S24**(-1)*S25 - 64.D0*S14**(-1)*S25**(-1)*S34 + 32.D
     &    0*S14**(-1) + 32.D0*S14**(-1)*S23*S25**(-1) + 32.D0*S14**(-1)
     &    *S15*S23**(-1)*S25**(-1)*S34 - 32.D0*S14**(-1)*S15*S23**(-1)
     &     + 128.D0*S14**(-1)*S15*S24**(-1)*S25**(-1)*S34 - 64.D0*
     &    S14**(-1)*S15*S24**(-1) + 48.D0*S14**(-1)*S15*S25**(-1) - 32.D
     &    0*S14**(-1)*S15**2*S24**(-1)*S25**(-1) - 32.D0*S15**(-1)*
     &    S23**(-1)*S24**(-1)*S34**2 - 64.D0*S15**(-1)*S23**(-1)*
     &    S24**(-1)*S25*S34 + 32.D0*S15**(-1)*S23**(-1)*S25**(-1)*
     &    S34**2 - 96.D0*S15**(-1)*S23**(-1)*S34 - 32.D0*S15**(-1)*
     &    S23**(-1)*S25 + 32.D0*S15**(-1)*S24**(-1)*S25**(-1)*S34**2 -
     &    96.D0*S15**(-1)*S24**(-1)*S34 - 32.D0*S15**(-1)*S24**(-1)*S25
     &     + 32.D0*S23**(-1)*S24**(-1)*S34 + 64.D0*S23**(-1)*S24**(-1)*
     &    S25 - 32.D0*S23**(-1)*S25**(-1)*S34 + 32.D0*S23**(-1) - 32.D0
     &    *S24**(-1)*S25**(-1)*S34 + 32.D0*S24**(-1) - 32.D0*S14*
     &    S15**(-1)*S23**(-1)*S25**(-1)*S34 )
      msq1 = msq1 + Pi**2*e**4*as*Cf/N * ( 32.D0*S14*
     &    S15**(-1)*S23**(-1) + 32.D0*S14*S23**(-1)*S25**(-1) - 32.D0*
     &    S13*S15**(-1)*S24**(-1)*S25**(-1)*S34 + 32.D0*S13*S15**(-1)*
     &    S24**(-1) + 32.D0*S13*S24**(-1)*S25**(-1) )
       msq(1) = msq1 


        msq1 =  Pi**2*e**4*as*Cf/N * (  - 64.D0*S13**(-1)*S14**(-1)*
     &    S23**(-1)*S24**(-1)*S34**3 - 64.D0*S13**(-1)*S14**(-1)*
     &    S23**(-1)*S34**2 - 64.D0*S13**(-1)*S14**(-1)*S24**(-1)*S34**2
     &     + 32.D0*S13**(-1)*S14**(-1)*S15*S23**(-1)*S25**(-1)*S34**2
     &     - 64.D0*S13**(-1)*S14**(-1)*S15*S23**(-1)*S34 + 32.D0*
     &    S13**(-1)*S14**(-1)*S15*S24**(-1)*S25**(-1)*S34**2 - 64.D0*
     &    S13**(-1)*S14**(-1)*S15*S24**(-1)*S34 + 16.D0*S13**(-1)*
     &    S14**(-1)*S15*S25**(-1)*S34 - 16.D0*S13**(-1)*S14**(-1)*S15
     &     - 32.D0*S13**(-1)*S14**(-1)*S15**2*S23**(-1) - 32.D0*
     &    S13**(-1)*S14**(-1)*S15**2*S24**(-1) - 16.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S24**(-1)*S34**3 - 16.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S24**(-1)*S25*S34**2 + 64.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S25**(-1)*S34**3 - 144.D0*S13**(-1)*
     &    S15**(-1)*S23**(-1)*S34**2 + 80.D0*S13**(-1)*S15**(-1)*
     &    S23**(-1)*S25*S34 - 32.D0*S13**(-1)*S15**(-1)*S23**(-1)*
     &    S25**2 )
     &  + Pi**2*e**4*as*Cf/N * (  - 16.D0*
     &    S13**(-1)*S15**(-1)*S24**(-1)*S34**2 - 16.D0*S13**(-1)*
     &    S15**(-1)*S24**(-1)*S25*S34 + 32.D0*S13**(-1)*S15**(-1)*
     &    S25**(-1)*S34**2 - 32.D0*S13**(-1)*S15**(-1)*S34 - 32.D0*
     &    S13**(-1)*S15**(-1)*S24*S25**(-1)*S34 + 32.D0*S13**(-1)*
     &    S15**(-1)*S24 - 48.D0*S13**(-1)*S23**(-1)*S24**(-1)*S34**2 -
     &    16.D0*S13**(-1)*S23**(-1)*S24**(-1)*S25*S34 - 32.D0*S13**(-1)
     &    *S23**(-1)*S24**(-1)*S25**2 - 128.D0*S13**(-1)*S23**(-1)*
     &    S25**(-1)*S34**2 + 64.D0*S13**(-1)*S23**(-1)*S34 - 64.D0*
     &    S13**(-1)*S23**(-1)*S25 - 64.D0*S13**(-1)*S24**(-1)*S34 - 32.D
     &    0*S13**(-1)*S24**(-1)*S25 - 64.D0*S13**(-1)*S25**(-1)*S34 +
     &    32.D0*S13**(-1) + 32.D0*S13**(-1)*S24*S25**(-1) + 128.D0*
     &    S13**(-1)*S15*S23**(-1)*S25**(-1)*S34 - 64.D0*S13**(-1)*S15*
     &    S23**(-1) + 32.D0*S13**(-1)*S15*S24**(-1)*S25**(-1)*S34 - 32.D
     &    0*S13**(-1)*S15*S24**(-1) + 48.D0*S13**(-1)*S15*S25**(-1) -
     &    32.D0*S13**(-1)*S15**2*S23**(-1)*S25**(-1) )
     &  + Pi**2*e**4*as*Cf/N * (  - 16.D0*
     &    S14**(-1)*S15**(-1)*S23**(-1)*S24**(-1)*S34**3 - 16.D0*
     &    S14**(-1)*S15**(-1)*S23**(-1)*S24**(-1)*S25*S34**2 - 16.D0*
     &    S14**(-1)*S15**(-1)*S23**(-1)*S34**2 - 16.D0*S14**(-1)*
     &    S15**(-1)*S23**(-1)*S25*S34 + 64.D0*S14**(-1)*S15**(-1)*
     &    S24**(-1)*S25**(-1)*S34**3 - 144.D0*S14**(-1)*S15**(-1)*
     &    S24**(-1)*S34**2 + 80.D0*S14**(-1)*S15**(-1)*S24**(-1)*S25*
     &    S34 - 32.D0*S14**(-1)*S15**(-1)*S24**(-1)*S25**2 + 32.D0*
     &    S14**(-1)*S15**(-1)*S25**(-1)*S34**2 - 32.D0*S14**(-1)*
     &    S15**(-1)*S34 - 32.D0*S14**(-1)*S15**(-1)*S23*S25**(-1)*S34
     &     + 32.D0*S14**(-1)*S15**(-1)*S23 - 48.D0*S14**(-1)*S23**(-1)*
     &    S24**(-1)*S34**2 - 16.D0*S14**(-1)*S23**(-1)*S24**(-1)*S25*
     &    S34 - 32.D0*S14**(-1)*S23**(-1)*S24**(-1)*S25**2 - 64.D0*
     &    S14**(-1)*S23**(-1)*S34 - 32.D0*S14**(-1)*S23**(-1)*S25 - 128.
     &    D0*S14**(-1)*S24**(-1)*S25**(-1)*S34**2 + 64.D0*S14**(-1)*
     &    S24**(-1)*S34 )
     &  + Pi**2*e**4*as*Cf/N * (  - 64.D0*
     &    S14**(-1)*S24**(-1)*S25 - 64.D0*S14**(-1)*S25**(-1)*S34 + 32.D
     &    0*S14**(-1) + 32.D0*S14**(-1)*S23*S25**(-1) + 32.D0*S14**(-1)
     &    *S15*S23**(-1)*S25**(-1)*S34 - 32.D0*S14**(-1)*S15*S23**(-1)
     &     + 128.D0*S14**(-1)*S15*S24**(-1)*S25**(-1)*S34 - 64.D0*
     &    S14**(-1)*S15*S24**(-1) + 48.D0*S14**(-1)*S15*S25**(-1) - 32.D
     &    0*S14**(-1)*S15**2*S24**(-1)*S25**(-1) - 32.D0*S15**(-1)*
     &    S23**(-1)*S24**(-1)*S34**2 - 64.D0*S15**(-1)*S23**(-1)*
     &    S24**(-1)*S25*S34 + 32.D0*S15**(-1)*S23**(-1)*S25**(-1)*
     &    S34**2 - 96.D0*S15**(-1)*S23**(-1)*S34 - 32.D0*S15**(-1)*
     &    S23**(-1)*S25 + 32.D0*S15**(-1)*S24**(-1)*S25**(-1)*S34**2 -
     &    96.D0*S15**(-1)*S24**(-1)*S34 - 32.D0*S15**(-1)*S24**(-1)*S25
     &     + 32.D0*S23**(-1)*S24**(-1)*S34 + 64.D0*S23**(-1)*S24**(-1)*
     &    S25 - 32.D0*S23**(-1)*S25**(-1)*S34 + 32.D0*S23**(-1) - 32.D0
     &    *S24**(-1)*S25**(-1)*S34 + 32.D0*S24**(-1) - 32.D0*S14*
     &    S15**(-1)*S23**(-1)*S25**(-1)*S34 )
     &  + Pi**2*e**4*as*Cf/N * ( 32.D0*S14*
     &    S15**(-1)*S23**(-1) + 32.D0*S14*S23**(-1)*S25**(-1) - 32.D0*
     &    S13*S15**(-1)*S24**(-1)*S25**(-1)*S34 + 32.D0*S13*S15**(-1)*
     &    S24**(-1) + 32.D0*S13*S24**(-1)*S25**(-1) )
       msq(2) = msq1 






c           print*,msq1
c
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
c         msq(2) = gs**2*qe**4*qu**2*((16*s14*s23)/(9.*s15*s34**2) + 
c     -    (16*s13*s24)/(9.*s15*s34**2) + (32*s23*s24)/(9.*s12*s34**2) - 
c     -    (16*s14*s23*s25)/(9.*s12*s15*s34**2) - 
c     -    (16*s13*s24*s25)/(9.*s12*s15*s34**2) + 
c     -    (16*s14*s35)/(9.*s12*s34**2) + (16*s24*s35)/(9.*s12*s34**2) - 
c     -    (16*s24*s35)/(9.*s15*s34**2) + 
c     -    (16*s14*s25*s35)/(9.*s12*s15*s34**2) + 
c     -    (32*s24*s25*s35)/(9.*s12*s15*s34**2) + 
c     -    (16*s13*s45)/(9.*s12*s34**2) + (16*s23*s45)/(9.*s12*s34**2) - 
c     -    (16*s23*s45)/(9.*s15*s34**2) + 
c     -    (16*s13*s25*s45)/(9.*s12*s15*s34**2) + 
c     -    (32*s23*s25*s45)/(9.*s12*s15*s34**2) - 
c     -    (32*s35*s45)/(9.*s15*s34**2))
c
c         msq(3)= gs**2*qe**4*qu**2*((32*s13*s14)/(9.*s12*s34**2) +
c     -    (16*s14*s23)/(9.*s25*s34**2) -
c     -    (16*s14*s15*s23)/(9.*s12*s25*s34**2) +
c     -    (16*s13*s24)/(9.*s25*s34**2) -
c     -    (16*s13*s15*s24)/(9.*s12*s25*s34**2) +
c     -    (16*s14*s35)/(9.*s12*s34**2) + (16*s24*s35)/(9.*s12*s34**2) -
c     -    (16*s14*s35)/(9.*s25*s34**2) +
c     -    (32*s14*s15*s35)/(9.*s12*s25*s34**2) +
c     -    (16*s15*s24*s35)/(9.*s12*s25*s34**2) +
c     -    (16*s13*s45)/(9.*s12*s34**2) + (16*s23*s45)/(9.*s12*s34**2) -
c     -    (16*s13*s45)/(9.*s25*s34**2) +
c     -    (32*s13*s15*s45)/(9.*s12*s25*s34**2) +
c     -    (16*s15*s23*s45)/(9.*s12*s25*s34**2) -
c     -    (32*s35*s45)/(9.*s25*s34**2))
c
c            ! Total contribution from qg ang qq initiated process
c                print*," "
c       msq(1) = msq1   ! + msq2 !this msq2 is the contribution from qg initiaied process
cc      msq = msq1
c           msq(1) =     (16*as*Cf*ge**4*Pi**2*
c     -    (2*S13**2*S14*S23*(S15 + S25 - S34) - 
c     -      S15*(2*S15**2*(S23 + S24)*S25 + 
c     -         2*S15*S24*(2*S25 - S34)*S34 + 
c     -         4*S25*S34**2*(S23 + S24 + S34) + 
c     -         S15*S23*(S24*S25 - S24*S34 + 4*S25*S34 - 2*S34**2)) - 
c     -      S13*(2*S15**3*S23 - 2*S23**2*S24*S25 + 2*S23*S25**3 + 
c     -         S15**2*(S23*(-3*S24 + 4*S25 - 8*S34) + 
c     -            2*S24*(S25 - S34)) + 2*S23**2*S24*S34 + 
c     -         2*S23*S24*S25*S34 - 5*S23*S25**2*S34 + S24*S25**2*S34 - 
c     -         2*S23*S24*S34**2 + 9*S23*S25*S34**2 + S24*S25*S34**2 + 
c     -         S25**2*S34**2 - 4*S23*S34**3 + S25*S34**3 + 
c     -         2*S14**2*S24*(-S25 + S34) + 
c     -         2*S14*(S25*S34*(2*S25 + S34) + 
c     -            S23*(S25**2 + 3*S25*S34 - S34**2) + 
c     -            S24*(S25**2 + 3*S25*S34 - S34**2)) + 
c     -         S15*(-2*S14**2*S24 - 2*S23**2*S24 + 
c     -            S23*(-2*S24*S25 + 4*S25**2 + 4*S24*S34 - 4*S25*S34 + 
c     -               8*S34**2) - 
c     -            2*S14*(S23*(S25 - S34) + S24*(S25 - S34) + 
c     -               S25*(2*S25 + S34)) + 
c     -            S25*(2*S25**2 + S25*S34 + 3*S34**2 + 
c     -               2*S24*(S25 + 2*S34)))) - 
c     -      S14*(2*S15**3*S24 + 2*S24*S25**3 + 
c     -         S15**2*(4*S24*(S25 - 2*S34) + 
c     -            S23*(-3*S24 + 2*S25 - 2*S34)) - 5*S24*S25**2*S34 + 
c     -         9*S24*S25*S34**2 + S25**2*S34**2 - 4*S24*S34**3 + 
c     -         S25*S34**3 + S23*
c     -          (-2*S24**2*(S25 - S34) + 2*S24*(S25 - S34)*S34 + 
c     -            S25*S34*(S25 + S34)) + 
c     -         S15*(4*S24*(S25**2 - S25*S34 + 2*S34**2) + 
c     -            S25*(2*S25**2 + S25*S34 + 3*S34**2) - 
c     -            2*S23*(S24**2 + S24*(S25 - 2*S34) - S25*(S25 + 2*S34))
c     -            ))))/(N*S13*S14*S15*S23*S24*S25)
c

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



         msq(3)=64d0*gs**2*e**4*(2*s12**3*(s13 + s15 
     -           - s35)*(s23 + s25 - s35)*
     -       (2*(s13 + s15)*(s23 + s25) - (s13 + s15 + s23 + s25)*s35 +
     -         s35**2) - 2*s12**2*(s13 + s15 - s35)*(s23 + s25 - s35)*
     -       (s15*(s15 + s23 + s25)*(3*s23 + 2*s25) +
     -         s13**2*(2*s23 + 3*s25 - s35) -
     -         (s15**2 + 4*s15*(s23 + s25) + (s23 + s25)**2)*s35 +
     -         (s15 + s23 + s25)*s35**2 +
     -         s13*((s23 + s25)*(5*s15 + 2*s23 + 3*s25) -
     -            2*(s15 + 2*(s23 + s25))*s35 + s35**2)) -
     -      (s13 + s15 + s23 + s25 - s35)*
     -       (s13**4*s25*(s23 + s25 - s35) +
     -         s13**3*(s25*(2*s25 - 3*s35)*(s23 + s25 - s35) -
     -            s15*(s23**2 - s23*s25 + s25*s35)) +
     -         s13*(-(s25*(s25 - s35)**2*(s23 + s25 - s35)*s35) +
     -            s15**3*((s23 - s25)*s25 - s23*s35) +
     -            s15*(s23**4 + s23**3*s25 + s23*s25**3 + s25**4 -
     -               (s23 + s25)*(3*s23**2 - s23*s25 + 3*s25**2)*s35 +
     -               2*(s23**2 + s23*s25 + s25**2)*s35**2 -
     -               (s23 + s25)*s35**3) +
     -            2*s15**2*(s23**3 + s23**2*(s25 - s35) + s25**2*s35 +
     -               s23*s35*(-s25 + s35))) +
     -         s15*s23*(s23 + s25 - s35)*
     -          (s15**3 + s15**2*(2*s23 - 3*s35) - (s23 - s35)**2*s35 +
     -            s15*(s23**2 - s25**2 + s25*s35 + 3*s35**2 -
     -               s23*(s25 + 4*s35))) +
     -         s13**2*(-2*s15**2*(s23**2 + s25**2) +
     -            s25*(s23 + s25 - s35)*
     -             (-s23**2 + (s25 - 3*s35)*(s25 - s35) +
     -               s23*(-s25 + s35)) +
     -            2*s15*(s23*s25*(s25 - s35) + s23**2*s35 +
     -               s25*(s25**2 - s25*s35 + s35**2)))) +
     -      s12*(s13**4*(s23 + s25 - s35)*(2*s23 + 4*s25 - s35) +
     -         (s23 + s25 - s35)*
     -          (2*s15**2*(s15**2*(2*s23 + s25) +
     -               s15*s23*(3*s23 + 2*s25) +
     -               (s23 + s25)*(2*s23**2 + s23*s25 + s25**2)) -
     -            s15*(s15**3 + s15**2*(11*s23 + 5*s25) +
     -               (s23 + s25)*(5*s23**2 + 5*s23*s25 + 3*s25**2) +
     -               s15*(13*s23**2 + 11*s23*s25 + 3*s25**2))*s35 +
     -            (2*s15**3 + (s23 + s25)**3 + 5*s15**2*(2*s23 + s25) +
     -               2*s15*(4*s23**2 + 5*s23*s25 + 2*s25**2))*s35**2 -
     -            (s15**2 + (s23 + s25)**2 + s15*(3*s23 + 2*s25))*s35**3
     -            ) + s13**3*
     -          ((s23 + s25 - s35)*
     -             (4*s23*s25 + 6*s25**2 - 5*s23*s35 - 11*s25*s35 +
     -               2*s35**2) +
     -            s15*(6*s23**2 + 18*s23*s25 + 10*s25**2 - 11*s23*s35 -
     -               15*s25*s35 + 4*s35**2)) +
     -         s13**2*(2*s15**2*
     -             (5*s23**2 + 12*s23*s25 + 5*s25**2 -
     -               9*(s23 + s25)*s35 + 3*s35**2) +
     -            2*s15*(2*s23**3 + 8*s23**2*s25 + 11*s23*s25**2 +
     -               5*s25**3 -
     -               (9*s23**2 + 25*s23*s25 + 15*s25**2)*s35 +
     -               (11*s23 + 14*s25)*s35**2 - 3*s35**3) +
     -            (s23 + s25 - s35)*
     -             (2*s23**3 + s23**2*(4*s25 - 3*s35) +
     -               s23*(6*s25 - 5*s35)*(s25 - s35) +
     -               (s25 - s35)*(4*s25**2 - 9*s25*s35 + s35**2))) +
     -         s13*(s15**3*(10*s23**2 + 18*s23*s25 + 6*s25**2 -
     -               15*s23*s35 - 11*s25*s35 + 4*s35**2) +
     -            s15*(6*(s23 + s25)**2*(s23**2 + s23*s25 + s25**2) -
     -               10*(s23 + s25)*(2*s23**2 + 3*s23*s25 + 2*s25**2)*
     -                s35 + 2*(13*s23**2 + 25*s23*s25 + 13*s25**2)*
     -                s35**2 - 15*(s23 + s25)*s35**3 + 2*s35**4) -
     -            (s23 + s25 - s35)*s35*
     -             (3*s23**3 + s23**2*(8*s25 - 4*s35) +
     -               s25*(5*s25 - 3*s35)*(s25 - s35) +
     -               2*s23*(5*s25**2 - 5*s25*s35 + s35**2)) +
     -            2*s15**2*(5*s23**3 + s23**2*(11*s25 - 15*s35) +
     -               (2*s25 - 3*s35)*(s25**2 - 3*s25*s35 + s35**2) +
     -               s23*(8*s25**2 - 25*s25*s35 + 14*s35**2)))))/
     -  (s13*s15*s23*s25*(s13 + s15 - s35)**2*(s23 + s25 - s35)**2)
         msq(3)=msq(3)/72d0
c           print*,msq(1)
c           print*,msq(2)
c           print*,msq(3)
c           print*," "
      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c     [u U -> a a]  Born 
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
      s13 =  -2.0d0*dot(p1,p3) ! t
c      s23 =  -2.0d0*dot(p1,p4) ! u
      s12 =   2.0d0*dot(p1,p2) ! s
      s23 =  -(s13+s12) ! u
      qu2 = 1d0!4d0/9d0
      xn = 3d0

      Born_uU2eE= CF*(e**4*(s13/s23 + s23/s13)/xn)

       return
       end
c---------------------------------------------------------------------
