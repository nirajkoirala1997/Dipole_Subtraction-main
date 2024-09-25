c---------------------------------------------------------------------
c       This file contains      
c     1.  Born / reduced_Born gg2aa in BSM 
c     2.  1real DrellYan ( qq ,qg ) channels 
c     3.  1real gg2aa  contribution to be filled
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     [g g -> a a]  Born BSM
c--------------------------------------------------------------------o
       function Born_gg2aa(k,p1,p2,p3,p4)
         implicit double precision (a-h,o-z)
         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
         parameter(PI=3.141592653589793238D0)
         common/usedalpha/AL,ge
c         common/energy/s12

         e= DSQRT(ge*4.d0*PI)

         s12 = 2d0*dot(p1,p2)
         rp34  = dsqrt(s12)

         model = 1    ! ADD model
         call coupfact(model,rp34,AK2D,AK2DINTF)

         s =  2d0*dot(p1,p2)
         t = -2d0*dot(p1,p3)
         u = -2d0*dot(p1,p4)
c         u = -(s+t)
c Reduced born used different kinematics, Colour factor multiply here.
         IF(k .eq. 0)  CF =  1d0               !Leading Order k=0
         IF(k .eq. 1)  CF = -4d0/3d0               !leg 1 reduced born k=1 
         IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2 reduced born k=2

c         IF(k .eq. 1)  CF = -1d0               !leg 1 reduced born k=1 
c         IF(k .eq. 2)  CF = -1d0               !Leg 2 reduced born k=2

         Born_gg2aa= CF*(u**4 + t**4)*AK2D**2/16.d0/8d0
c	print*,"Born1:",Born_gg2aa
c         Born_gg2aa= (4 (t**4 + u**4))/(t + u + xmg)^2
c	print*,"Born2:",Born_gg2aa
c	print*," "

       return
       end
c---------------------------------------------------------------------

       subroutine amp_mat_r(p,msq)
       implicit double precision (a-h,o-z)
       dimension p(0:3,1:5),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
     . ,a(1:100)
       parameter(Pi=3.141592653589793238D0)
       double precision msq(5),msq1,msq2,lambda
       common/usedalpha/AL,ge
       common/scales/xinvmass

       call p2dtop1d_5(p,p1,p2,p3,p4,p5)

c       s12=2.d0*dot(p1,p2)
c       s13=2.d0*dot(p1,p3)
c       s14=2.d0*dot(p1,p4)
c       s15=2.d0*dot(p1,p5)
c       s23=2.d0*dot(p2,p3)
c       s24=2.d0*dot(p2,p4)
c       s25=2.d0*dot(p2,p5)
c       s34=2.d0*dot(p3,p4)
c       s35=2.d0*dot(p3,p5)
c       s45=2.d0*dot(p4,p5)

       P12=2d0*dot(p1,p2)
       P13=-2d0*dot(p1,p3)
       P14=-2d0*dot(p1,p4)
       P23=-2d0*dot(p2,p3)
       P24=-2d0*dot(p2,p4)

       P34=2d0*dot(p3,p4)
       P35=2d0*dot(p3,p5)
       P45=2d0*dot(p4,p5)
       P15=-2d0*dot(p1,p5)
       P25=-2d0*dot(p2,p5)



c       P35=P12+P14+P24
c       P45=P12+P13+P23
c       P34=P12-P35-P45
c       P15=P23+P24+P34
c       P25=P13+P14+P34


       qe=DSQRT(ge*4.d0*PI)
c      gs=DSQRT(AL*4.d0*PI)
       qu =1d0! 2d0/3d0

       aem = ge
       model = 1
       lambda = 0.226d0 
       rp34  = dsqrt(P34)
       call coupfact(model,rp34,AK2D,AK2DINTF)

c________________________________________________________________________c 
       N=3
       Cf=(N*N-1.0D0)/2.D0/N
       Tf=1.0D0/2.0D0
       N2m1=N*N-1.0D0
       e=DSQRT(4.0D0*PI*aem)

       nf = 5
       xnf=nf
       xmur = xinvmass

       as = AL/4.0D0/PI


c       print*,"P12=",P12   
c       print*,"P13=",P13  
c       print*,"P14=",P14  
c       print*,"P23=",P23  
c       print*,"P24=",P24  
c       print*,"P34=",P34  
c       print*,"P35=",P35  
c       print*,"P45=",P45  
c       print*,"P15=",P15  
c       print*,"P25=",P25  
c
c       print*,"   N =",N
c       print*,"N2m1 =",N2m1
c       print*,"  as =",as
c       print*,"ak2D =",ak2D


c      SGRgg =
c     &   N*Pi**2*as/N2m1*ak2D**2 * ( 4.D0*P12**(-2)*P34**5 - 4.D0
c     &    *P12**(-2)*P25*P34**4 + 8.D0*P12**(-2)*P24*P25*P34**3 - 8.D0
c     &    *P12**(-2)*P24*P25**2*P34**2 + 4.D0*P12**(-2)*P24**2*P34**3
c     &     - 4.D0*P12**(-2)*P24**2*P25*P34**2 - 8.D0*P12**(-2)*P23**2*
c     &    P25*P34**2 + 12.D0*P12**(-2)*P23**2*P25**2*P34 - 4.D0*
c     &    P12**(-2)*P23**2*P25**3 - 12.D0*P12**(-2)*P15*P34**4 + 8.D0*
c     &    P12**(-2)*P15*P25*P34**3 - 24.D0*P12**(-2)*P15*P24*P25*P34**2
c     &     + 16.D0*P12**(-2)*P15*P24*P25**2*P34 - 12.D0*P12**(-2)*P15*
c     &    P24**2*P34**2 + 8.D0*P12**(-2)*P15*P24**2*P25*P34 + 16.D0*
c     &    P12**(-2)*P15*P23*P34**3 - 24.D0*P12**(-2)*P15*P23*P25*P34**2
c     &     + 8.D0*P12**(-2)*P15*P23*P25**2*P34 + 8.D0*P12**(-2)*P15*
c     &    P23**2*P34**2 - 4.D0*P12**(-2)*P15*P23**2*P25**2 + 12.D0*
c     &    P12**(-2)*P15**2*P34**3 - 4.D0*P12**(-2)*P15**2*P25*P34**2 +
c     &    24.D0*P12**(-2)*P15**2*P24*P25*P34 - 8.D0*P12**(-2)*P15**2*
c     &    P24*P25**2 + 12.D0*P12**(-2)*P15**2*P24**2*P34 - 4.D0*
c     &    P12**(-2)*P15**2*P24**2*P25 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 40.D0*
c     &    P12**(-2)*P15**2*P23*P34**2 + 40.D0*P12**(-2)*P15**2*P23*P25*
c     &    P34 - 8.D0*P12**(-2)*P15**2*P23*P25**2 - 12.D0*P12**(-2)*
c     &    P15**2*P23**2*P34 + 4.D0*P12**(-2)*P15**2*P23**2*P25 - 4.D0*
c     &    P12**(-2)*P15**3*P34**2 - 8.D0*P12**(-2)*P15**3*P24*P25 - 4.D0
c     &    *P12**(-2)*P15**3*P24**2 + 32.D0*P12**(-2)*P15**3*P23*P34 -
c     &    16.D0*P12**(-2)*P15**3*P23*P25 + 4.D0*P12**(-2)*P15**3*P23**2
c     &     - 8.D0*P12**(-2)*P15**4*P23 + 8.D0*P12**(-2)*P14*P34**4 - 8.D
c     &    0*P12**(-2)*P14*P25*P34**3 + 16.D0*P12**(-2)*P14*P23*P34**3
c     &     - 24.D0*P12**(-2)*P14*P23*P25*P34**2 + 8.D0*P12**(-2)*P14*
c     &    P23*P25**2*P34 - 24.D0*P12**(-2)*P14*P15*P34**3 + 16.D0*
c     &    P12**(-2)*P14*P15*P25*P34**2 - 40.D0*P12**(-2)*P14*P15*P23*
c     &    P34**2 + 40.D0*P12**(-2)*P14*P15*P23*P25*P34 - 8.D0*P12**(-2)
c     &    *P14*P15*P23*P25**2 + 24.D0*P12**(-2)*P14*P15**2*P34**2 - 8.D0
c     &    *P12**(-2)*P14*P15**2*P25*P34 + 32.D0*P12**(-2)*P14*P15**2*
c     &    P23*P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 16.D0*
c     &    P12**(-2)*P14*P15**2*P23*P25 - 8.D0*P12**(-2)*P14*P15**3*P34
c     &     - 8.D0*P12**(-2)*P14*P15**3*P23 + 4.D0*P12**(-2)*P14**2*
c     &    P34**3 - 4.D0*P12**(-2)*P14**2*P25*P34**2 - 12.D0*P12**(-2)*
c     &    P14**2*P15*P34**2 + 8.D0*P12**(-2)*P14**2*P15*P25*P34 + 12.D0
c     &    *P12**(-2)*P14**2*P15**2*P34 - 4.D0*P12**(-2)*P14**2*P15**2*
c     &    P25 - 4.D0*P12**(-2)*P14**2*P15**3 - 8.D0*P12**(-2)*P13*P24*
c     &    P34**3 + 8.D0*P12**(-2)*P13*P24*P25*P34**2 + 24.D0*P12**(-2)*
c     &    P13*P15*P24*P34**2 - 16.D0*P12**(-2)*P13*P15*P24*P25*P34 - 24.
c     &    D0*P12**(-2)*P13*P15**2*P24*P34 + 8.D0*P12**(-2)*P13*P15**2*
c     &    P24*P25 + 8.D0*P12**(-2)*P13*P15**3*P24 + 2.D0*P12**(-1)*
c     &    P15**(-1)*P34**5 + 6.D0*P12**(-1)*P15**(-1)*P25**2*P34**3 -
c     &    16.D0*P12**(-1)*P15**(-1)*P24*P34**4 - 12.D0*P12**(-1)*
c     &    P15**(-1)*P24*P25*P34**3 + 36.D0*P12**(-1)*P15**(-1)*P24*
c     &    P25**2*P34**2 + 16.D0*P12**(-1)*P15**(-1)*P24*P25**3*P34 + 32.
c     &    D0*P12**(-1)*P15**(-1)*P24**2*P34**3 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 100.D0*
c     &    P12**(-1)*P15**(-1)*P24**2*P25*P34**2 + 104.D0*P12**(-1)*
c     &    P15**(-1)*P24**2*P25**2*P34 - 24.D0*P12**(-1)*P15**(-1)*
c     &    P24**3*P34**2 + 48.D0*P12**(-1)*P15**(-1)*P24**3*P25*P34 + 8.D
c     &    0*P12**(-1)*P15**(-1)*P24**4*P34 + 40.D0*P12**(-1)*P15**(-1)*
c     &    P23*P34**4 - 8.D0*P12**(-1)*P15**(-1)*P23*P25*P34**3 - 56.D0*
c     &    P12**(-1)*P15**(-1)*P23*P25**2*P34**2 + 24.D0*P12**(-1)*
c     &    P15**(-1)*P23*P24*P34**3 - 72.D0*P12**(-1)*P15**(-1)*P23*P24*
c     &    P25*P34**2 - 16.D0*P12**(-1)*P15**(-1)*P23*P24**2*P34**2 + 96.
c     &    D0*P12**(-1)*P15**(-1)*P23**2*P34**3 - 72.D0*P12**(-1)*
c     &    P15**(-1)*P23**2*P25*P34**2 - 60.D0*P12**(-1)*P15**(-1)*
c     &    P23**2*P25**2*P34 + 24.D0*P12**(-1)*P15**(-1)*P23**2*P24*
c     &    P34**2 - 72.D0*P12**(-1)*P15**(-1)*P23**2*P24*P25*P34 - 16.D0
c     &    *P12**(-1)*P15**(-1)*P23**2*P24**2*P34 + 72.D0*P12**(-1)*
c     &    P15**(-1)*P23**3*P34**2 - 64.D0*P12**(-1)*P15**(-1)*P23**3*
c     &    P25*P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * ( 16.D0*
c     &    P12**(-1)*P15**(-1)*P23**4*P34 + 2.D0*P12**(-1)*P25**(-1)*
c     &    P34**5 - 16.D0*P12**(-1)*P34**4 + 16.D0*P12**(-1)*P25*P34**3
c     &     - 2.D0*P12**(-1)*P25**2*P34**2 - 16.D0*P12**(-1)*P24*
c     &    P25**(-1)*P34**4 - 16.D0*P12**(-1)*P24*P34**3 + 28.D0*
c     &    P12**(-1)*P24*P25*P34**2 + 32.D0*P12**(-1)*P24*P25**2*P34 +
c     &    32.D0*P12**(-1)*P24**2*P25**(-1)*P34**3 - 124.D0*P12**(-1)*
c     &    P24**2*P34**2 + 136.D0*P12**(-1)*P24**2*P25*P34 - 24.D0*
c     &    P12**(-1)*P24**3*P25**(-1)*P34**2 + 56.D0*P12**(-1)*P24**3*
c     &    P34 + 8.D0*P12**(-1)*P24**4*P25**(-1)*P34 + 40.D0*P12**(-1)*
c     &    P23*P25**(-1)*P34**4 - 48.D0*P12**(-1)*P23*P34**3 - 72.D0*
c     &    P12**(-1)*P23*P25*P34**2 + 80.D0*P12**(-1)*P23*P25**2*P34 +
c     &    24.D0*P12**(-1)*P23*P24*P25**(-1)*P34**3 - 96.D0*P12**(-1)*
c     &    P23*P24*P34**2 + 72.D0*P12**(-1)*P23*P24*P25*P34 - 16.D0*
c     &    P12**(-1)*P23*P24**2*P25**(-1)*P34**2 + 16.D0*P12**(-1)*P23*
c     &    P24**2*P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * ( 96.D0*
c     &    P12**(-1)*P23**2*P25**(-1)*P34**3 - 168.D0*P12**(-1)*P23**2*
c     &    P34**2 + 20.D0*P12**(-1)*P23**2*P25*P34 - 4.D0*P12**(-1)*
c     &    P23**2*P25**2 + 24.D0*P12**(-1)*P23**2*P24*P25**(-1)*P34**2
c     &     - 72.D0*P12**(-1)*P23**2*P24*P34 - 16.D0*P12**(-1)*P23**2*
c     &    P24**2*P25**(-1)*P34 + 72.D0*P12**(-1)*P23**3*P25**(-1)*
c     &    P34**2 - 104.D0*P12**(-1)*P23**3*P34 + 16.D0*P12**(-1)*P23**4
c     &    *P25**(-1)*P34 + 28.D0*P12**(-1)*P15*P34**3 - 20.D0*P12**(-1)
c     &    *P15*P24*P25**(-1)*P34**3 + 44.D0*P12**(-1)*P15*P24*P34**2 +
c     &    24.D0*P12**(-1)*P15*P24*P25*P34 - 4.D0*P12**(-1)*P15*P24**2*
c     &    P25**(-1)*P34**2 + 36.D0*P12**(-1)*P15*P24**2*P34 + 8.D0*
c     &    P12**(-1)*P15*P24**3*P25**(-1)*P34 - 56.D0*P12**(-1)*P15*P23*
c     &    P25**(-1)*P34**3 - 8.D0*P12**(-1)*P15*P23*P34**2 + 120.D0*
c     &    P12**(-1)*P15*P23*P25*P34 - 8.D0*P12**(-1)*P15*P23*P25**2 -
c     &    24.D0*P12**(-1)*P15*P23*P24*P25**(-1)*P34**2 + 72.D0*
c     &    P12**(-1)*P15*P23*P24*P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * ( 16.D0*
c     &    P12**(-1)*P15*P23*P24**2*P25**(-1)*P34 - 96.D0*P12**(-1)*P15*
c     &    P23**2*P25**(-1)*P34**2 + 100.D0*P12**(-1)*P15*P23**2*P34 -
c     &    40.D0*P12**(-1)*P15*P23**3*P25**(-1)*P34 + 2.D0*P12**(-1)*
c     &    P15**2*P25**(-1)*P34**3 + 2.D0*P12**(-1)*P15**2*P34**2 - 8.D0
c     &    *P12**(-1)*P15**2*P25*P34 + 4.D0*P12**(-1)*P15**2*P24*
c     &    P25**(-1)*P34**2 - 16.D0*P12**(-1)*P15**2*P24*P25 - 4.D0*
c     &    P12**(-1)*P15**2*P24**2*P25**(-1)*P34 - 4.D0*P12**(-1)*P15**2
c     &    *P24**2 + 40.D0*P12**(-1)*P15**2*P23*P25**(-1)*P34**2 + 40.D0
c     &    *P12**(-1)*P15**2*P23*P34 - 24.D0*P12**(-1)*P15**2*P23*P25 +
c     &    36.D0*P12**(-1)*P15**2*P23**2*P25**(-1)*P34 + 4.D0*P12**(-1)*
c     &    P15**2*P23**2 - 8.D0*P12**(-1)*P15**3*P34 - 8.D0*P12**(-1)*
c     &    P15**3*P24 - 16.D0*P12**(-1)*P15**3*P23*P25**(-1)*P34 - 16.D0
c     &    *P12**(-1)*P15**3*P23 + 12.D0*P12**(-1)*P14*P15**(-1)*P25*
c     &    P34**3 + 12.D0*P12**(-1)*P14*P15**(-1)*P25**2*P34**2 - 24.D0*
c     &    P12**(-1)*P14*P15**(-1)*P24*P34**3 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 24.D0*
c     &    P12**(-1)*P14*P15**(-1)*P24*P25*P34**2 + 48.D0*P12**(-1)*P14*
c     &    P15**(-1)*P24*P25**2*P34 + 56.D0*P12**(-1)*P14*P15**(-1)*P23*
c     &    P34**3 - 16.D0*P12**(-1)*P14*P15**(-1)*P23*P25*P34**2 + 8.D0*
c     &    P12**(-1)*P14*P15**(-1)*P23*P25**2*P34 + 72.D0*P12**(-1)*P14*
c     &    P15**(-1)*P23**2*P34**2 - 24.D0*P12**(-1)*P14*P15**(-1)*
c     &    P23**2*P25*P34 + 16.D0*P12**(-1)*P14*P15**(-1)*P23**3*P34 + 8.
c     &    D0*P12**(-1)*P14*P34**3 + 20.D0*P12**(-1)*P14*P25*P34**2 - 24.
c     &    D0*P12**(-1)*P14*P24*P25**(-1)*P34**3 - 48.D0*P12**(-1)*P14*
c     &    P24*P34**2 + 72.D0*P12**(-1)*P14*P24*P25*P34 + 56.D0*
c     &    P12**(-1)*P14*P23*P25**(-1)*P34**3 - 96.D0*P12**(-1)*P14*P23*
c     &    P34**2 + 80.D0*P12**(-1)*P14*P23*P25*P34 + 72.D0*P12**(-1)*
c     &    P14*P23**2*P25**(-1)*P34**2 - 24.D0*P12**(-1)*P14*P23**2*P34
c     &     + 16.D0*P12**(-1)*P14*P23**3*P25**(-1)*P34 + 20.D0*P12**(-1)
c     &    *P14*P15*P25**(-1)*P34**3 + 28.D0*P12**(-1)*P14*P15*P34**2 -
c     &    24.D0*P12**(-1)*P14*P15*P24*P25**(-1)*P34**2 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * ( 24.D0*
c     &    P12**(-1)*P14*P15*P24*P34 - 16.D0*P12**(-1)*P14*P15*P23*
c     &    P25**(-1)*P34**2 + 72.D0*P12**(-1)*P14*P15*P23*P34 - 8.D0*
c     &    P12**(-1)*P14*P15*P23*P25 - 4.D0*P12**(-1)*P14*P15**2*
c     &    P25**(-1)*P34**2 - 8.D0*P12**(-1)*P14*P15**2*P25 - 16.D0*
c     &    P12**(-1)*P14*P15**2*P23*P25**(-1)*P34 - 8.D0*P12**(-1)*P14*
c     &    P15**2*P23 - 8.D0*P12**(-1)*P14*P15**3 + 8.D0*P12**(-1)*
c     &    P14**2*P15**(-1)*P34**3 + 20.D0*P12**(-1)*P14**2*P15**(-1)*
c     &    P25*P34**2 + 8.D0*P12**(-1)*P14**2*P15**(-1)*P25**2*P34 - 32.D
c     &    0*P12**(-1)*P14**2*P15**(-1)*P24*P34**2 + 32.D0*P12**(-1)*
c     &    P14**2*P15**(-1)*P24*P25*P34 + 24.D0*P12**(-1)*P14**2*
c     &    P15**(-1)*P23*P34**2 + 24.D0*P12**(-1)*P14**2*P15**(-1)*P23*
c     &    P25*P34 + 24.D0*P12**(-1)*P14**2*P15**(-1)*P23**2*P34 + 8.D0*
c     &    P12**(-1)*P14**2*P25**(-1)*P34**3 + 20.D0*P12**(-1)*P14**2*
c     &    P34**2 + 16.D0*P12**(-1)*P14**2*P25*P34 - 32.D0*P12**(-1)*
c     &    P14**2*P24*P25**(-1)*P34**2 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * ( 32.D0*
c     &    P12**(-1)*P14**2*P24*P34 + 24.D0*P12**(-1)*P14**2*P23*
c     &    P25**(-1)*P34**2 + 24.D0*P12**(-1)*P14**2*P23*P34 + 24.D0*
c     &    P12**(-1)*P14**2*P23**2*P25**(-1)*P34 + 20.D0*P12**(-1)*
c     &    P14**2*P15*P25**(-1)*P34**2 + 12.D0*P12**(-1)*P14**2*P15*P34
c     &     - 4.D0*P12**(-1)*P14**2*P15**2*P25**(-1)*P34 - 4.D0*
c     &    P12**(-1)*P14**2*P15**2 + 8.D0*P12**(-1)*P14**3*P15**(-1)*
c     &    P34**2 + 16.D0*P12**(-1)*P14**3*P15**(-1)*P25*P34 + 16.D0*
c     &    P12**(-1)*P14**3*P15**(-1)*P23*P34 + 8.D0*P12**(-1)*P14**3*
c     &    P25**(-1)*P34**2 + 24.D0*P12**(-1)*P14**3*P34 + 16.D0*
c     &    P12**(-1)*P14**3*P23*P25**(-1)*P34 + 8.D0*P12**(-1)*P14**3*
c     &    P15*P25**(-1)*P34 + 8.D0*P12**(-1)*P14**4*P15**(-1)*P34 + 8.D0
c     &    *P12**(-1)*P14**4*P25**(-1)*P34 + 8.D0*P12**(-1)*P13*
c     &    P15**(-1)*P24*P34**3 + 80.D0*P12**(-1)*P13*P15**(-1)*P24*P25*
c     &    P34**2 - 64.D0*P12**(-1)*P13*P15**(-1)*P24*P25**2*P34 + 120.D0
c     &    *P12**(-1)*P13*P15**(-1)*P24**2*P34**2 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 48.D0*
c     &    P12**(-1)*P13*P15**(-1)*P24**2*P25*P34 + 16.D0*P12**(-1)*P13*
c     &    P15**(-1)*P24**3*P34 - 48.D0*P12**(-1)*P13*P15**(-1)*P23*
c     &    P34**3 - 48.D0*P12**(-1)*P13*P15**(-1)*P23*P24*P34**2 - 48.D0
c     &    *P12**(-1)*P13*P15**(-1)*P23**2*P34**2 - 48.D0*P12**(-1)*P13*
c     &    P15**(-1)*P23**2*P24*P34 + 8.D0*P12**(-1)*P13*P24*P25**(-1)*
c     &    P34**3 + 104.D0*P12**(-1)*P13*P24*P34**2 - 200.D0*P12**(-1)*
c     &    P13*P24*P25*P34 + 120.D0*P12**(-1)*P13*P24**2*P25**(-1)*
c     &    P34**2 - 120.D0*P12**(-1)*P13*P24**2*P34 + 16.D0*P12**(-1)*
c     &    P13*P24**3*P25**(-1)*P34 - 48.D0*P12**(-1)*P13*P23*P25**(-1)*
c     &    P34**3 + 96.D0*P12**(-1)*P13*P23*P34**2 - 48.D0*P12**(-1)*P13
c     &    *P23*P24*P25**(-1)*P34**2 + 48.D0*P12**(-1)*P13*P23*P24*P34
c     &     - 48.D0*P12**(-1)*P13*P23**2*P25**(-1)*P34**2 + 48.D0*
c     &    P12**(-1)*P13*P23**2*P34 - 48.D0*P12**(-1)*P13*P23**2*P24*
c     &    P25**(-1)*P34 - 16.D0*P12**(-1)*P13*P15*P24*P25**(-1)*P34**2
c     &     - 144.D0*P12**(-1)*P13*P15*P24*P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 72.D0*
c     &    P12**(-1)*P13*P15*P24**2*P25**(-1)*P34 + 96.D0*P12**(-1)*P13*
c     &    P15*P23*P25**(-1)*P34**2 - 48.D0*P12**(-1)*P13*P15*P23*P34 +
c     &    48.D0*P12**(-1)*P13*P15*P23*P24*P25**(-1)*P34 + 48.D0*
c     &    P12**(-1)*P13*P15*P23**2*P25**(-1)*P34 + 8.D0*P12**(-1)*P13*
c     &    P15**2*P24*P25**(-1)*P34 + 8.D0*P12**(-1)*P13*P15**2*P24 - 48.
c     &    D0*P12**(-1)*P13*P15**2*P23*P25**(-1)*P34 + 32.D0*P12**(-1)*
c     &    P13*P14*P15**(-1)*P24*P34**2 - 32.D0*P12**(-1)*P13*P14*
c     &    P15**(-1)*P24*P25*P34 + 32.D0*P12**(-1)*P13*P14*P24*P25**(-1)
c     &    *P34**2 - 32.D0*P12**(-1)*P13*P14*P24*P34 - 24.D0*P12**(-1)*
c     &    P13**2*P15**(-1)*P24*P34**2 + 48.D0*P12**(-1)*P13**2*
c     &    P15**(-1)*P24*P25*P34 - 24.D0*P12**(-1)*P13**2*P24*P25**(-1)*
c     &    P34**2 + 120.D0*P12**(-1)*P13**2*P24*P34 + 72.D0*P12**(-1)*
c     &    P13**2*P15*P24*P25**(-1)*P34 + 32.D0*P12**(-1)*P13**2*P14*
c     &    P15**(-1)*P24*P34 + 32.D0*P12**(-1)*P13**2*P14*P24*P25**(-1)*
c     &    P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * ( 40.D0*
c     &    P15**(-1)*P25**(-1)*P34**5 - 34.D0*P15**(-1)*P34**4 + 26.D0*
c     &    P15**(-1)*P25*P34**3 - 4.D0*P15**(-1)*P25**2*P34**2 - 80.D0*
c     &    P15**(-1)*P24*P25**(-1)*P34**4 + 128.D0*P15**(-1)*P24*P34**3
c     &     + 20.D0*P15**(-1)*P24*P25*P34**2 + 16.D0*P15**(-1)*P24*
c     &    P25**2*P34 + 12.D0*P15**(-1)*P24**2*P25**(-1)*P34**3 - 20.D0*
c     &    P15**(-1)*P24**2*P34**2 + 104.D0*P15**(-1)*P24**2*P25*P34 + 8.
c     &    D0*P15**(-1)*P24**3*P25**(-1)*P34**2 + 48.D0*P15**(-1)*P24**3
c     &    *P34 + 8.D0*P15**(-1)*P24**4*P25**(-1)*P34 + 104.D0*P15**(-1)
c     &    *P23*P25**(-1)*P34**4 - 128.D0*P15**(-1)*P23*P34**3 - 48.D0*
c     &    P15**(-1)*P23*P25*P34**2 - 24.D0*P15**(-1)*P23*P24*P25**(-1)*
c     &    P34**3 - 72.D0*P15**(-1)*P23*P24*P34**2 - 16.D0*P15**(-1)*P23
c     &    *P24**2*P25**(-1)*P34**2 + 96.D0*P15**(-1)*P23**2*P25**(-1)*
c     &    P34**3 - 168.D0*P15**(-1)*P23**2*P34**2 - 60.D0*P15**(-1)*
c     &    P23**2*P25*P34 - 24.D0*P15**(-1)*P23**2*P24*P25**(-1)*P34**2
c     &     - 72.D0*P15**(-1)*P23**2*P24*P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 16.D0*
c     &    P15**(-1)*P23**2*P24**2*P25**(-1)*P34 + 40.D0*P15**(-1)*
c     &    P23**3*P25**(-1)*P34**2 - 64.D0*P15**(-1)*P23**3*P34 + 16.D0*
c     &    P15**(-1)*P23**4*P25**(-1)*P34 - 50.D0*P25**(-1)*P34**4 + 48.D
c     &    0*P34**3 - 14.D0*P25*P34**2 - 16.D0*P24*P25**(-1)*P34**3 - 16.
c     &    D0*P24*P34**2 + 16.D0*P24*P25*P34 - 44.D0*P24**2*P25**(-1)*
c     &    P34**2 + 32.D0*P24**2*P34 + 8.D0*P24**3*P25**(-1)*P34 - 56.D0
c     &    *P23*P25**(-1)*P34**3 + 104.D0*P23*P34**2 + 80.D0*P23*P25*P34
c     &     + 24.D0*P23*P24*P25**(-1)*P34**2 + 72.D0*P23*P24*P34 + 16.D0
c     &    *P23*P24**2*P25**(-1)*P34 - 24.D0*P23**2*P25**(-1)*P34**2 +
c     &    72.D0*P23**2*P34 - 40.D0*P23**3*P25**(-1)*P34 + 34.D0*P15*
c     &    P25**(-1)*P34**3 - 10.D0*P15*P34**2 + 28.D0*P15*P24*P25**(-1)
c     &    *P34**2 - 8.D0*P15*P24*P34 - 4.D0*P15*P24**2*P25**(-1)*P34 +
c     &    24.D0*P15*P23*P34 - 8.D0*P15*P23*P25 + 36.D0*P15*P23**2*
c     &    P25**(-1)*P34 - 8.D0*P15**2*P25**(-1)*P34**2 - 8.D0*P15**2*
c     &    P34 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 8.D0*P15**2
c     &    *P24 - 16.D0*P15**2*P23*P25**(-1)*P34 - 8.D0*P15**2*P23 + 88.D
c     &    0*P14*P15**(-1)*P25**(-1)*P34**4 - 32.D0*P14*P15**(-1)*P34**3
c     &     + 28.D0*P14*P15**(-1)*P25*P34**2 - 120.D0*P14*P15**(-1)*P24*
c     &    P25**(-1)*P34**3 + 72.D0*P14*P15**(-1)*P24*P34**2 + 48.D0*P14
c     &    *P15**(-1)*P24*P25*P34 + 168.D0*P14*P15**(-1)*P23*P25**(-1)*
c     &    P34**3 - 32.D0*P14*P15**(-1)*P23*P34**2 + 8.D0*P14*P15**(-1)*
c     &    P23*P25*P34 + 120.D0*P14*P15**(-1)*P23**2*P25**(-1)*P34**2 -
c     &    24.D0*P14*P15**(-1)*P23**2*P34 + 16.D0*P14*P15**(-1)*P23**3*
c     &    P25**(-1)*P34 - 56.D0*P14*P25**(-1)*P34**3 + 40.D0*P14*P34**2
c     &     - 24.D0*P14*P24*P25**(-1)*P34**2 + 24.D0*P14*P24*P34 - 80.D0
c     &    *P14*P23*P25**(-1)*P34**2 + 64.D0*P14*P23*P34 + 20.D0*P14*P15
c     &    *P25**(-1)*P34**2 - 16.D0*P14*P15*P23*P25**(-1)*P34 - 8.D0*
c     &    P14*P15**2 + 84.D0*P14**2*P15**(-1)*P25**(-1)*P34**3 + 4.D0*
c     &    P14**2*P15**(-1)*P34**2 + 8.D0*P14**2*P15**(-1)*P25*P34 - 32.D
c     &    0*P14**2*P15**(-1)*P24*P25**(-1)*P34**2 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * ( 32.D0*P14**2*
c     &    P15**(-1)*P24*P34 + 72.D0*P14**2*P15**(-1)*P23*P25**(-1)*
c     &    P34**2 + 24.D0*P14**2*P15**(-1)*P23*P34 + 24.D0*P14**2*
c     &    P15**(-1)*P23**2*P25**(-1)*P34 - 20.D0*P14**2*P25**(-1)*
c     &    P34**2 + 8.D0*P14**2*P34 - 4.D0*P14**2*P15*P25**(-1)*P34 + 40.
c     &    D0*P14**3*P15**(-1)*P25**(-1)*P34**2 + 16.D0*P14**3*P15**(-1)
c     &    *P34 + 16.D0*P14**3*P15**(-1)*P23*P25**(-1)*P34 + 8.D0*P14**3
c     &    *P25**(-1)*P34 + 8.D0*P14**4*P15**(-1)*P25**(-1)*P34 + 48.D0*
c     &    P13*P15**(-1)*P24*P25**(-1)*P34**3 + 16.D0*P13*P15**(-1)*P24*
c     &    P34**2 - 64.D0*P13*P15**(-1)*P24*P25*P34 + 120.D0*P13*
c     &    P15**(-1)*P24**2*P25**(-1)*P34**2 - 48.D0*P13*P15**(-1)*
c     &    P24**2*P34 + 16.D0*P13*P15**(-1)*P24**3*P25**(-1)*P34 - 48.D0
c     &    *P13*P15**(-1)*P23*P25**(-1)*P34**3 - 48.D0*P13*P15**(-1)*P23
c     &    *P24*P25**(-1)*P34**2 - 48.D0*P13*P15**(-1)*P23**2*P25**(-1)*
c     &    P34**2 - 48.D0*P13*P15**(-1)*P23**2*P24*P25**(-1)*P34 - 32.D0
c     &    *P13*P24*P25**(-1)*P34**2 )
c      SGRgg = SGRgg + N*Pi**2*as/N2m1*ak2D**2 * (  - 136.D0*P13*
c     &    P24*P34 - 72.D0*P13*P24**2*P25**(-1)*P34 + 96.D0*P13*P23*
c     &    P25**(-1)*P34**2 + 48.D0*P13*P23*P24*P25**(-1)*P34 + 48.D0*
c     &    P13*P23**2*P25**(-1)*P34 + 8.D0*P13*P15*P24*P25**(-1)*P34 -
c     &    48.D0*P13*P15*P23*P25**(-1)*P34 + 32.D0*P13*P14*P15**(-1)*P24
c     &    *P25**(-1)*P34**2 - 32.D0*P13*P14*P15**(-1)*P24*P34 + 72.D0*
c     &    P13**2*P15**(-1)*P24*P25**(-1)*P34**2 + 48.D0*P13**2*
c     &    P15**(-1)*P24*P34 + 72.D0*P13**2*P24*P25**(-1)*P34 + 32.D0*
c     &    P13**2*P14*P15**(-1)*P24*P25**(-1)*P34 )
c      msq(4) = SGRgg

c      print*,msq(4)
c      stop
c________________________________________________________________________c 
c
c
        aN2m1=8d0
         aNn =3d0
cc         
cc
         do i=1,99
           a(i) = 0d0
         enddo
cc
      include "amplitude.h"
      msq(4) = SGRgg1
cc      print*,"mod:",SGRgg1
cc
cc         do i=1,99
c           print*,i,a(i)
c         enddo
c      stop

      return
      end
c---------------------------------------------------------------------

