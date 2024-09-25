      function flo2_Vir(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),coef(2,-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      parameter (cf=4d0/3d0)
      parameter (zeta2=1.64493406684823d0)
      parameter (N=3)
      parameter (Nf=1)
      parameter (Ca=3)
      parameter (Tf=0.5D0) 
      parameter (EulerGamma = 0.5772156649015328606065120d0)

      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      external Born_gg2aa
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        scale  = xinvmass
        if ( scale .ge. xlow .and. scale .le. xhigh) then 
             
              xmuf=scale
              xmur=scale
              xmu2=xmuf**2

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)

              e= DSQRT(ge*4.d0*PI)
              gs=DSQRT(Al*4.d0*PI)
              qu2=1d0                           ! charge is considered in lum.f

              s12 =  2d0*dot(p1,p2)
               t  = -2d0*dot(p1,p3)
               u  = -2d0*dot(p1,p4)

             Born = Born_gg2aa(0,p1,p2,p3,p4)

c**********[ From Slicing Code ] **********o
        vir = Born*( -N*203.d0/9.d0
     &            +1.*Tf*70.d0/9.d0
     &            +N*zeta2*8.0d0 )
c**********[ From Slicing Ends ] **********o
        model =1
        rp34 = dsqrt(s12)
        call coupfact(model,rp34,AK2D,AK2DINTF)
c	print*,AK2D,EulerGamma,Pi,s,t,u
c	stop

          eikonal_born = -0.00005425347222222222d0*
     &  (AK2D**2*AL*(CA*(264d0*t**3*u + 462d0*t**2*u**2 + 264d0*t*u**3 +
     &  t**4*(-235d0 + 132d0*EulerGamma + 12d0*Pi**2 - 264d0*Log(2d0)-
     &       132d0*Log(Pi)) + 
     &    u**4*(-235d0 + 132d0*EulerGamma + 12d0*Pi**2 -264d0*Log(2d0)-
     &       132d0*Log(Pi))) - 
     & 4d0*Nf*TR*(24d0*t**3*u + 42d0*t**2*u**2 + 24d0*t*u**3 + 
     &    t**4*(-17d0 + 12d0*EulerGamma - 24d0*Log(2d0)-12d0*Log(Pi)) +
     &    u**4*(-17d0 + 12d0*EulerGamma - 24d0*Log(2d0)-12d0*Log(Pi)))-
     & 6d0*CF*(6d0*t**2*u**2*
     &     (6d0 + 7d0*EulerGamma - 14d0*Log(2d0) - 7d0*Log(Pi)) + 
     &    6d0*t**3*u*(3d0 + 4d0*EulerGamma - 8d0*Log(2d0)- 4d0*Log(Pi))+
     &    6*t*u**3*(3d0 + 4d0*EulerGamma - 8d0*Log(2d0) - 4d0*Log(Pi))+
     &    t**4*(3d0 + 6d0*EulerGamma**2 - 5d0*Pi**2 + 24d0*Log(2d0)**2-
     &       15d0*Log(Pi) + 6d0*Log(Pi)**2 + 
     &       6d0*Log(2d0)*(-5d0 + 4d0*Log(Pi)) - 
     &       3d0*EulerGamma*(-5d0 + Log(256d0) + 4d0*Log(Pi))) + 
     &    u**4*(3d0 + 6d0*EulerGamma**2 - 5d0*Pi**2+24d0*Log(2d0)**2 - 
     &       15d0*Log(Pi) + 6d0*Log(Pi)**2 + 
     &       6d0*Log(2d0)*(-5d0 + 4d0*Log(Pi)) - 
     &       3d0*EulerGamma*(-5d0 + Log(256d0) + 4d0*Log(Pi)))) + 
     & 6d0*(-22d0*CA*(t**4 + u**4) + 8d0*Nf*TR*(t**4 + u**4) + 
     &    3d0*CF*(8d0*t**3*u + 14d0*t**2*u**2 + 8d0*t*u**3d0 + 
     &       t**4*(5d0 + 4d0*EulerGamma - 8d0*Log(2d0)- 4d0*Log(Pi)) + 
     &       u**4*(5d0 + 4d0*EulerGamma - 8d0*Log(2d0)- 4d0*Log(Pi))))*
     &  Log(xmu2/s12) - 36d0*CF*(t**4 + u**4)*Log(xmu2/s12)**2))/Pi 



        sig= xl(4)*(Vir + eikonal_born)  
c       sig= xl(1)*SumI(1)   

         xnorm=hbarc2/16d0/pi/(xa*xb*s)
         wgt=xnorm*sig*vwgt
         flo2_Vir=wgt/vwgt/2d0/eps

            return
       else                  
        flo2_Vir=0d0
        return
       endif
c       else
c        flo2_Vir=0d0
c       endif
      end

c---------------------------------------------------------------------

C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer n4,ndim,ncomp,nvec,core,iter,userdata 
      real*8 xx(ndim) ,f(ncomp),weight,flo2_Vir
      common/countc/n4
      external flo2_Vir
      f(1) = flo2_Vir(xx,weight)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
