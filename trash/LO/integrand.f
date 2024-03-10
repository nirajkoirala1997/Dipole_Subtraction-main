C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer n4,ndim,ncomp,nvec,core,iter,userdata 
      real*8 xx(ndim) ,f(ncomp),weight,flo2_Vir
      common/countc/n4
      real *8 flo2_LO
      external flo2_LO
      f(1) = flo2_LO(xx,weight)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo2_LO(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
c      parameter (hbarc2=0.3894d9)
      parameter (hbarc2=389.3856741d6)
      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
c      common/renor_scale/scale
      common/usedalpha/AL,ge
      external Born_uU2eE
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
        eps = 0.5d0
       xlow = xq - eps
      xhigh = xq + eps

c      write(*,*)'Q,Qmin,Qmax =', xq,xlow,xhigh
c       write(*,*)'Energy = ', dsqrt(s)

       call kinvar2(yy,xinvmass,p1,p2,p3,p4)

       scale2 = 2.0d0*dot(p3,p4)
       scale1 = dsqrt(scale2)

        if ( scale1 .ge. xlow .and. scale1 .le. xhigh) then 
             
c              xmuf=scale
c              xmur=scale
c              xmu2=xmuf**2
c              xmuf=xq
c              xmur=xq

c              call pdf(xa,xmuf,f1)
c              call pdf(xb,xmuf,f2)
c              call setlum(f1,f2,xl)
c              AL = alphasPDF(xmur)
c                ALSWZ=0.120d0
c                XMT = 172.5d0
c                call InitAlphaS(1, 1.0D0, 91.1876D0, ALSWZ,
c     &                  1.4D0, 4.75D0, XMT )
c              AL = alphaS(xmur)

c              sig= xl(1)*Born_uU2eE(0,p1,p2,p3,p4)
c               sig = 100.0d0

c              xnorm=hbarc2/16d0/pi/(xa*xb*s)
c              wgt=xnorm*sig*vwgt
c              flo2_LO=wgt/vwgt/2d0/eps
              flo2_LO=yy(1)*yy(2)*yy(3)
       else                  
        flo2_LO=0d0
       endif
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
       e= DSQRT(ge*4.d0*PI)
       gs=DSQRT(Al*4.d0*PI)

      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0 
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2    
      s13 =  2.0d0*dot(p1,p3) ! t
      s23 =  2.0d0*dot(p2,p3) ! u
      s12 =  2.0d0*dot(p1,p2) ! s
      XNC = 1/4d0
      xnorm =1d0
      qu2 = 1d0!4d0/9d0
c       2 added because of u u~ symm
      Born_uU2eE= 1d0 * (2d0*e**4*qu2*(-2*s13*s23 + s12*(s13 +
     .            s23)))/(3d0*s12**2)
       return
       end
c---------------------------------------------------------------------

