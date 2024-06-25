C Cuba Specific function.      
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
C -------------------------------------------------------------------- C        
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(2),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .             , sig(1:5),SumD(1:2)
      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4,unphy
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/amass/am1,am2,am3,am4,am5
      common/energy/s
      common/set/set1
      common/distribution/xq
      common/countc/n4
      common/usedalpha/AL,ge
      common/scales/xinvmass
      common/momenta5/p1,p2,p3,p4,p5
      external dipole_gq_q




      xa = xx(1)
      xb = xx(2)
      rsp = dsqrt(xa*xb*s)

      xmz = 91.1876d0

      ipass1 = 0

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

        fnlo3 = 0


c       call kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)


c------------------------------------------------------------------------
      call kinvar3_slicing(xx,xxjac,xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2, !
     &     r34,r35,r45,Et5,QT)                                          ! 
c------------------------------------------------------------------------
c      call cuts3(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,                    !
c     &     r34,r35,r45,Et5,QT,ipass)                                    ! 
c------------------------------------------------------------------------

c        stop






c        call kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)
        call cuts3(p1,p2,p3,p4,p5,icuts)

        if (unphy .eq. 0 .and. icuts .eq. 1) then ! with zero unphysical PS points proceed
     
        scale = xinvmass
        if ( scale .ge. xlow .and. scale .le. xhigh) then

          xmuf=scale
          xmur=scale

          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call setlum(f1,f2,xl)

c          AL = alphasPDF(xmur)

          call p1dtop2d_5(p1,p2,p3,p4,p5,p)
          call mat_amp_r(p,sig)
c
c	print*,"xx(1)=",xx(1)
c	print*,"xx(2)=",xx(2)
c	print*,"xx(3)=",xx(3)
c	print*,"xx(4)=",xx(4)
c	print*,"xx(5)=",xx(5)
c	print*,"xx(6)=",xx(6)

c	print*," " 

c        print*,"P1:",p1
c        print*,"P2:",p2
c        print*,"P3:",p3
c        print*,"P4:",p4
c        print*,"P5:",p5


          SumD(1) = dipole_gg_g(1,p) + dipole_gg_g(2,p)

c          sigma = xl(4)*( sig(1)-SumD(1) )
          sigma = xl(4)*sig(4)
c          sigma = xl(4)*SumD(1)
          pi_1 = 0.5d0*rsp
          flux = 4d0*pi_1*rsp
          xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
          wgt=xxjac*xnorm*sigma*weight
          fnlo3=wgt/weight/2d0/eps

            xnorm=hbarc2*xx(3)/64/4/4/(pi**4)
            wgt=xxjac*xnorm*sigma*weight
            fnlo3=wgt/weight/2.d0/eps




c	print*,"Integrand:",fnlo3,xl(4)
c	stop
         endif
       endif

	if ( fnlo3 .ne. fnlo3) then
	 fnlo3 = 0d0
	endif

151   return
      end
c---------------------------------------------------------------------
