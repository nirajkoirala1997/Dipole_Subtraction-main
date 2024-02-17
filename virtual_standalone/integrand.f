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
      dimension f1(-6:6),f2(-6:6),xl(1)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      character*50 name
      common/pdfname/name
      common/leg_choice/leg
      s = 1000**2

      call kinvar2(yy,xinvmass,p1,p2,p3,p4)
      call p1d_to_p2d_4(p1,p2,p3,p4,p)
      xa     = yy(1)
      xb     = yy(2)
      AL     = 0.118d0
      ge=0.007547169811320755d0
      Al=0.118d0
      e= DSQRT(ge*4.d0*PI)
      gs=DSQRT(Al*4.d0*PI)
      qu2=4d0/9d0


      scale  = xinvmass
      if( scale .gt. 100d0) then
            xmuf=scale
            xmur=scale
       call pdf(xa,xmuf,f1)
       call pdf(xb,xmuf,f2)
       call setlum(f1,f2,xl)

        s= 2d0*dot(p1,p2)
        t=-2d0*dot(p2,p4)
        u=-2d0*dot(p2,p3)

c        VIR=(4*Al*e**4*qu2*(5*(-6 + 5*Pi**2)*s**2 + 
c     -     2*(-48 + 25*Pi**2)*s*t + 2*(-48 + 25*Pi**2)*t**2 - 
c     -     6*(s**2 + 6*s*t + 6*t**2)*Log((xmuf**2/s)) - 
c     -     6*(s**2 + 2*s*t + 2*t**2)*Log((xmuf**2/s))**2))/(3.*Pi*s**2)

        VIR= (32*AL*(-2 + Pi**2)*e**4*qu2*(s**2 +
     .          2*s*t + 2*t**2))/(Pi*s**2)
        call Iterm(p,coef,SumI)

c        print*,SumI

        sig= xl(1)*(Vir+Born_uU2eE(0,p1,p2,p3,p4)*SumI(0))
c        sig= xl(1)*Born_uU2eE(0,p1,p2,p3,p4)*SumI(0)

         xnorm=hbarc2/16d0/pi/s
         wgt=xnorm*sig*vwgt
         flo2_Vir=wgt/vwgt

         if( flo2_Vir .ne. flo2_Vir) flo2_Vir=0
                                ! There is a factor of 1/2
         else                   ! from xeps interval
            flo2_Vir=0d0
         endif
         return
      end

c---------------------------------------------------------------------

c--------------------------------------------------------------------o
       function Born_uU2eE(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       ge=0.007547169811320755d0
       Al=0.118d0
       e= DSQRT(ge*4.d0*PI)
       gs=DSQRT(Al*4.d0*PI)
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2
      s13 =  2.0d0*dot(p1,p3) ! t
      s23 =  2.0d0*dot(p2,p3) ! u
      s12 =  2.0d0*dot(p1,p2) ! s
      qu2 = 4d0/9d0

      Born_uU2eE= CF*(2*e**4*qu2*(-2*s13*s23 + s12*(s13 +
     .            s23)))/(3d0*s12**2)
       return
       end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
