C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer ndim,ncomp,nvec,core,iter,userdata 
      real*8 xx(ndim) ,f(ncomp),weight,fnlo3
      external fnlo3
      f(1) = fnlo3(xx,weight)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c     Our vegas specific
      double precision function fnlo3(xx,weight)
      implicit double precision(a-h,o-z) 
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      dimension  xx(6),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(2),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      common/amass/am1,am2,am3,am4,am5
      common/energy/s

      call kinvar3(xx,xxjac,p1,p2,p3,p4,p5)
      s12=am1**2 + am2**2 + 2d0*dot(p1,p2)
      s35=am3**2 + am5**2 + 2d0*dot(p3,p5)
      s45=am4**2 + am5**2 + 2d0*dot(p4,p5)
      s34=dot(p3,p4)

      s23=2d0*dot(p2,p3)
      s13=2d0*dot(p1,p3)
      s14=2d0*dot(p1,p4)
      s15=2d0*dot(p1,p5)
      s25=2d0*dot(p2,p5)
      s24=2d0*dot(p2,p4)

     
      xa=xx(1)
      xb=xx(2)
      sp = xa*xb*s
      rsp= dsqrt(sp)
      scale = dsqrt(dabs(s34))!xinvmass


      ipass = 0

      if ( scale .ge. 140.0d0 ) ipass=1
       fnlo3 = 0.0d0
      if (p1(0) .eq. 0d0) goto 151
      if ( ipass .eq. 1 ) then
         xmuf=scale
         xmur=scale
        call pdf(xa,xmuf,f1)
        call pdf(xb,xmuf,f2)
        call setlum(f1,f2,xl)
        call p1dtop2d_5(p1,p2,p3,p4,p5,p)

        call  uu2ee_r(p,sig)
        SumD=dipole_uU_g(1,p)+dipole_uU_g(2,p)

        sig2 = xl(1)*(sig - SumD)


         pi_1 = 0.5d0*rsp
         flux = 4d0*pi_1*rsp
         xnorm=hbarc2/8d0/(2d0*Pi)**4/flux

         wgt=xxjac*xnorm*sig2*weight
         fnlo3=wgt/weight
c         if (fnlo3 .ne. fnlo3) fnlo3 = 0d0
        endif
151      return
      end

