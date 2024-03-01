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
      function flo2_Vir(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      external Born_uU2eE
       
      rs  = dsqrt(s)
      xa     = yy(1)
      xb     = yy(2)

      rsp = dsqrt(xa*xb*s)
        
      ipass = 0
        eps = 0.5d0*2d0
       xlow = xq - eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

      if (rsp .gt. xcut) then

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
        scale  = xinvmass
        if ( scale .ge. xlow .and. scale .le. xhigh) then 
             
              xmuf=scale
              xmur=scale
              xmu2=xmuf**2
c              xmuf=xq
c              xmur=xq

              call pdf(xa,xmuf,f1)
              call pdf(xb,xmuf,f2)
              call setlum(f1,f2,xl)
              AL = alphasPDF(xmur)

c              ge=1d0/128d0
              e= DSQRT(ge*4.d0*PI)
              gs=DSQRT(Al*4.d0*PI)
              qu2=4d0/9d0
              qu2=1d0
              Cf = 4d0/3d0

              s12= 2d0*dot(p1,p2)
              t=-2d0*dot(p2,p4)
              u=-2d0*dot(p2,p3)

c
              VIR= gs**2*8*e**4*qu2*(5*(-6 + 11*Pi**2)*s12**2 +
     -        2*(-48 + 55*Pi**2)*s12*t + 2*(-48 + 55*Pi**2)*t**2 -
     -        6*(s12**2 + 6*s12*t + 6*t**2)*Log(xmu2/s12) -
     -        6*(s12**2 + 2*s12*t + 2*t**2)*Log(xmu2/s12)**2)/
     -          (24*Pi**2*s12**2)
               VIR = VIR/36d0

               gammaq = 3d0/2d0*Cf
               xkq     = Cf*(7d0/2d0-pi*pi/6d0)
           EulerGamma = 0.5772156649015328606065120d0

          SumI(0) = 
     -      (-6*Al*CF*EulerGamma**2 - 12*Al*CF*gammaq + 
     -     12*Al*CF*EulerGamma*gammaq - 12*Al*CF*xkq + 5*Al*CF*Pi**2 + 
     -     12*Al*CF*EulerGamma*Log(4*Pi) - 12*Al*CF*gammaq*Log(4*Pi) - 
     -     6*Al*CF*Log(4*Pi)**2 + 12*Al*CF*EulerGamma*Log(xmu2/s12) - 
     -     12*Al*CF*gammaq*Log(xmu2/s12) - 
     -     12*Al*CF*Log(4*Pi)*Log(xmu2/s12) - 6*Al*CF*Log(xmu2/s12)**2)/
     -      (24.*Pi)
          SumI(0) = SumI(0)*Born_uU2eE(0,p1,p2,p3,p4)!/9d0 ! 9 for colour  avg
c              call p1d_to_p2d_4(p1,p2,p3,p4,p)
c             call Iterm(p,coef,SumI)
c


              sig= xl(1)*(Vir + SumI(0))

              xnorm=hbarc2/16d0/pi/(xa*xb*s)
              wgt=xnorm*sig*vwgt
              flo2_Vir=wgt/vwgt/2d0/eps
            return
       else                  
        flo2_Vir=0d0
        return
       endif
       else
        flo2_Vir=0d0
       endif
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
c       Al=0.118d0
       e= DSQRT(ge*4.d0*PI)
       gs=DSQRT(Al*4.d0*PI)
c       write(*,*)'e = ',e
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0 
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2    
      s13 =  2.0d0*dot(p1,p3) ! t
      s23 =  2.0d0*dot(p2,p3) ! u
      s12 =  2.0d0*dot(p1,p2) ! s
      XNC = 1/4d0
      xnorm =1d0
      qu2 = 1d0!4d0/9d0

      Born_uU2eE=2d0* CF*(2*e**4*qu2*(-2*s13*s23 + s12*(s13 +
     .            s23)))/(3d0*s12**2)
       return
       end
c---------------------------------------------------------------------

