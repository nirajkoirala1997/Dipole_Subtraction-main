c---------------------------------------------------------------------
c       This file contains      
c     1.  Born / reduced_Born gg2h in EFT(mass of top infinite approxomation) 
c     2.  
c     3.  
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     [g g -> h]  Born Effective field theory
c--------------------------------------------------------------------o
       function Born_gg2h(k,p1,p2,p3)
         implicit double precision (a-h,o-z)
         dimension p1(0:3),p2(0:3),p3(0:3)
         parameter(PI=3.141592653589793238D0)
         common/usedalpha/AL,ge
         common/param2/xmur
         common/mass/amh
c         common/energy/s12

         e= DSQRT(ge*4.d0*PI)

         s12 = 2d0*dot(p1,p2)
         rp34  = dsqrt(s12)

         
c  Reduced born used different kinematics, Colour factor multiply here.
         IF(k .eq. 0)  CF =  1d0               !Leading Order k=0
         IF(k .eq. 1)  CF = -4d0/3d0               !leg 1 reduced born k=1 
         IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2 reduced born k=2

c         IF(k .eq. 1)  CF = -1d0               !leg 1 reduced born k=1 
c         IF(k .eq. 2)  CF = -1d0               !Leg 2 reduced born k=2

          NA = 8
          AS = alphasPDF(xmur)
c          AS = AS/4d0/PI

            v = 246d0
           ch = -4d0*AS/3d0/v 
          ch2 = ch * ch
c          Born_gg2h= 16d0*AL**2*amh**4d0/(3d0*PI*v)**2/NA
c           Born_gg2h= AS**2/72d0/PI/v**2/8d0

c 3manually.frm :amp=  1/36*PI^-4*v^-2*s12^2*AL^2    .or.  amp= 1d0/2d0*s12**2*ch2
c          Born_gg2H =  1d0/2d0*s12**2*ch2
          Born_gg2H = CF * 0.5d0*amh**4*ch2/16d0

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


c       qe=DSQRT(ge*4.d0*PI)
cc      gs=DSQRT(AL*4.d0*PI)
c       qu =1d0! 2d0/3d0
c
c       aem = ge
c       model = 1
c       lambda = 0.226d0 
c       rp34  = dsqrt(P34)
c
c       N=3
c       Cf=(N*N-1.0D0)/2.D0/N
c       Tf=1.0D0/2.0D0
c       N2m1=N*N-1.0D0
c       e=DSQRT(4.0D0*PI*aem)
c
c       nf = 5
c       xnf=nf
c       xmur = xinvmass
c
c       as = AL/4.0D0/PI
c
cc	calculate real matrix amplitued g g > H

      return
      end
c---------------------------------------------------------------------

