c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function flo1_LO(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)      ! in Pb
c      parameter (hbarc2=0.3894d12)    ! in Fb
      common/leg_choice/leg
      common/energy/s
      common/distribution/xq
      common/renor_scale/scale
      common/usedalpha/AL,ge
      common/mass/amh
      common/param2/xmur
      external Born_gg2H
       
        xa = yy(1)
c        xb = yy(2) 
        xb = amh**2d0/xa/s
        sp = xa*xb*S
       rsp = dsqrt(sp)

      call kinvar1(xa,xb,p1,p2,p3)

      scale = 2d0 * dot(p1,p2) 
	eps = 0.05d0
	Q2 = dsqrt(scale)
	Qmin = amh - eps
	Qmax = amh + eps

c	if(Q2 .ge. Qmin .and. Q2 .le. Qmax) then
c      xmuf=scale / 2d0
      xmuf= 62.5d0
      xmur=scale / 2d0
      xmu2=xmuf**2

      call pdf(xa,xmuf,f1)
      call pdf(xb,xmuf,f2)
      call setlum(f1,f2,xl)

      sig= xl(2)*Born_gg2H(0,p1,p2,p3)

       pi_1 = 0.5d0*rsp
       flux = 4d0*pi_1*rsp
       xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
c       xnorm=hbarc2/8d0/(2d0*Pi)**4/flux/2d0/eps

        flo1_LO  = xnorm*sig
c      else
c        flo1_LO= 0d0
c      endif

      return
      end

c---------------------------------------------------------------------
cC ~~~~~~~~~~~~~~~~~~~C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,yy,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer ndim,ncomp,nvec,core,iter,userdata 
c      real*8 xx(ndim) ,f(ncomp),weight,fnlo3,yy(10)
      real*8 yy(2) ,f(ncomp),weight,flo1_LO,xx(10)
      external flo1_LO
       xx(1) = yy(1)
       xx(2) = yy(2)
      f(1)=flo1_LO(xx,weight)
c	print*,"Fun "
c	print*,xx
c      f= xx(1)*xx(2)*xx(3)*xx(4)*xx(5)*xx(6)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
