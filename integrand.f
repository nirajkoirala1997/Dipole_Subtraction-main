C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer n4,ndim,ncomp,nvec,core,iter,userdata 
      real*8 xx(ndim) ,f(ncomp),weight,fnlo3
      common/countc/n4
      external fnlo3
      f(1) = fnlo3(xx,weight)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c      Our vegas specific
      double precision function fnlo3(xx,weight)
      implicit double precision(a-h,o-z) 
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(2),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      common/amass/am1,am2,am3,am4,am5
      common/energy/s
      common/set/set1
      common/countc/n4

      xa = xx(1)
      xb = xx(2)
      rsp = dsqrt(xa*xb*s)

      ipass1 = 0

      xq = 200.0d0
      eps = 0.5d0
      xlow = xq -eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

      if (rsp .gt. xcut) then
      ipass1 = 1
      endif

      if (ipass1 .eq.1) then
      call kinvar3(xx,xxjac,p1,p2,p3,p4,p5)


c      if (p4(0) .le.0.0d0) then 
c      fnlo3 = 0.0d0
c      goto 151
c      return
c      endif

      s12=am1**2 + am2**2 + 2d0*dot(p1,p2)
      s35=am3**2 + am5**2 + 2d0*dot(p3,p5)
      s45=am4**2 + am5**2 + 2d0*dot(p4,p5)
      s34=2d0*dot(p3,p4)

      s23=2d0*dot(p2,p3)
      s13=2d0*dot(p1,p3)
      s14=2d0*dot(p1,p4)
      s15=2d0*dot(p1,p5)
      s25=2d0*dot(p2,p5)
      s24=2d0*dot(p2,p4)

     
      sp = xa*xb*s
      rsp= dsqrt(sp)

      scale = dsqrt(s34)!xinvmass

      ipass = 0
      fnlo3 = 0

      if ( scale .ge. xlow .and. scale .le. xhigh) then
      ipass=1
      endif


      if ( ipass .eq. 1 ) then

c      e4L = p4(0)
c      print*,e4L

c      xmuf=scale
c      xmur=scale

         xmuf=xq
         xmur=xq

        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)
        call p1dtop2d_5(p1,p2,p3,p4,p5,p)

        call  uu2ee_r(p,sig)
c        SumD=dipole_uU_g(1,p)+dipole_uU_g(2,p)
        dip1=dipole_uU_g(1,p)
        dip2=dipole_uU_g(2,p)
        SumD  = dip1 + dip2
c        print*,"Sig :",sig
c        if (dip1 .lt. -11110d0) print*,"SumD:",SumD,dip1,dip2,sig
c        if (dip2 .lt. -11110d0) print*,"SumD:",SumD,dip1,dip2,sig
c        print*," "

        sig2 = xl(1)*(sig - SumD)

         pi_1 = 0.5d0*rsp
         flux = 4d0*pi_1*rsp
         xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
         wgt=xxjac*xnorm*sig2*weight
         fnlo3=wgt/weight
        endif
        endif
151      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo2_PK(yy,vwgt)
      implicit double precision (a-h,o-z)
      integer leg
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(1)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/leg_choice/leg
      call kinvar2(yy,xinvmass,p1,p2,p3,p4)
      call p1d_to_p2d_4(p1,p2,p3,p4,p)
      xa     = yy(1)
      xb     = yy(2)
       x     = yy(4)
      scale  = xinvmass
      if( scale .gt. 100d0) then
            xmuf=scale
            xmur=scale
       call pdf(xa,xmuf,f1)
       call pdf(xb,xmuf,f2)
       call setlum(f1,f2,xl)

       call getPK(leg,x,xmuf,p,SumP,SumK)

        do i=0,3
        xp1(i) = x*p1(i)
        xp2(i) = x*p2(i)
        enddo
       if (leg .eq. 1) Born(leg) = Born_uU2eE_PK(0,xp1,p2,p3,p4)
       if (leg .eq. 2) Born(leg) = Born_uU2eE_PK(0,p1,xp2,p3,p4)

       PK =SumP+SumK

       sig = xl(1)*Born(leg)*PK

c         xnorm=hbarc2/16d0/pi/(s*xa*xb)
         xnorm=hbarc2/16d0/pi/s
         wgt=xnorm*sig*vwgt
         flo2_PK=wgt/vwgt

         return              ! There is a factor of 1/2
         else                   ! from xeps interval
            flo2_PK=0d0
           return
         endif
         return
      end

c---------------------------------------------------------------------
c     [u U -> e E]  Born 
c--------------------------------------------------------------------o
       function Born_uU2eE_PK(k,p1,p2,p3,p4)
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

      Born_uU2eE_PK= CF*(2*e**4*qu2*(-2*s13*s23 + s12*(s13 +
     .            s23)))/(3d0*s12**2)
       return
       end

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
        
        VIR=(4*Al*e**4*qu2*(5*(-6 + 5*Pi**2)*s**2 + 
     -     2*(-48 + 25*Pi**2)*s*t + 2*(-48 + 25*Pi**2)*t**2 - 
     -     6*(s**2 + 6*s*t + 6*t**2)*Log((xmuf**2/s)) - 
     -     6*(s**2 + 2*s*t + 2*t**2)*Log((xmuf**2/s))**2))/(3.*Pi*s**2)
        call Iterm(p,coef,SumI)

c        sig= xl(1)*Born_uU2eE(0,p1,p2,p3,p4)*(Vir+SumI(0))
        sig= xl(1)*(Vir+Born_uU2eE(0,p1,p2,p3,p4)*SumI(0))
        
c        sig=Born_uU2eE_PK(0,p1,p2,p3,p4)*SumI(0)
c        print*,Vir,sig,SumI(0)
c        print*,"  "
c        stop
         xnorm=hbarc2/16d0/pi/s
         wgt=xnorm*sig*vwgt
         flo2_Vir=wgt/vwgt
c         if( flo2_Vir .ne. flo2_Vir) flo2_Vir=0
         return              ! There is a factor of 1/2
         else                   ! from xeps interval
            flo2_Vir=0d0
            return
         endif
         return
      end

c---------------------------------------------------------------------
c
