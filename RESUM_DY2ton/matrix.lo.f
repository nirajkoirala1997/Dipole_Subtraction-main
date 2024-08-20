c-------------------------------------------c
c-     2 -> 2 matrix element squared       -c
c-------------------------------------------c      

c     SM AND BSM
      subroutine setsigma(s,t,u,sigma)

c        use handyG

!zz sm
c        implicit double precision (a)

C        implicit none
        integer i,nf
        integer ORD,counts,m11switch
        real *8 qm,zm,am,hm,wm
        real *8 N,pi,pi2,aS,asp,ALEM,ALEM2,GMZ,GMW
        real *8 ghqq,azqq,gaqq,ghaz,ghzz,gw,ctw,Gf,ca,sw2,st2,stw
        real *8 dwa,dwz,sigma,sqamph,vev,yb,yt,alphasPDF,xlambda
        real(KIND=16) xlogumz,xlogtmz
        real *8 s,t,u,xmz,xmz2,xmh,xmh2,mq(6),ghq(6),gphiq(6),azq(6),
     &  gphiz,gz2,cw2,xmw,xmw2,xmuf,xmur,xnf,alfas1,alfas2,xmb,xmt,xmc
        real *8 gau,gvu,gad,gvd
        real *8 x,y,z,b,sn,un,yy,zz
        real *8 kk,tlval,thval
        real *8 born,oneL,twoL11,twoL02a,twoL02b,twoL02c,twoL02d,m11
        real *8 w,Bf,cf,cf2,zeta3,a,bt0
        dimension sigma(10),w(3000),a(30)
        double complex ans(10)
c        common/nflavour/nf
c        common/mass/z,h
c        common/quarkmass/xmb,xmt
        COMMON/XXMH/XMUR,XMUF
        COMMON/ALS/AS
        COMMON/WENBERG/SW2,CW2,ALEM
        COMMON/XMASS/XMT,XMB,XMC
        COMMON/VMASS/XMZ,XMW,XMH,GMZ,GMW
        COMMON/NORDER/ORD
        COMMON/m11ONOFF/m11switch
        COMMON/Counter/counts
        DOUBLE PRECISION dilog
        double precision m0m0,m0m1,m1m1,m0m2a,m0m2b,m0m2c,m0m2d
        EXTERNAL dilog
        EXTERNAL cpp_wrapper

c        call clearcache

        xnf = 5d0
        ALEM2 = ALEM*ALEM
        Zeta3 = 1.2020569031595942853997381615114499907649862923405d0
        ctw = dsqrt(CW2)                            !Cos(thetaW)
        stw = dsqrt(sw2)

        gvu = 0.5d0*(0.5d0-2.0d0*sw2*2.0d0/3.0d0)/(stW*ctW)
        gvd = 0.5d0*(-0.5d0-2.0d0*sW2*(-1.0d0/3.0d0))/(stW*ctW)
        gau = 1.0d0/4d0/ctw/stw
        gad = gau

        PI = dacos(-1D0)
        pi2 = pi*pi
        un = u/xmz**2 
        sn = s/xmz**2
        
        b = dsqrt(1d0-4d0/sn)
        z = -un
        x = (1d0 - b)/(1d0 + b)
        y = (1d0 + x**2 - x*z)
        y = y/x
        yy = 1d0 + y
        zz = 1d0 + z
       !born factor
c        Bf = 4d0/9.0d0*ALEM2*Pi**2
        Bf = 3.0d0*1.0d0/4d0/9d0*(ALEM*4.0d0*Pi)**2.0d0
        cf = 4d0/3d0
        cf2 = cf*cf
        CA = 3d0
        BT0 = 11./3.D0 * CA - 2./3.D0 * xNF

        
        !Physical region
        !Kallen function
        kk = s*dsqrt( 1 - 4*xmz**2/s )
        tlval = 1d0/2d0*(2d0*xmz**2-s-kk)
        thval = 1d0/2d0*(2d0*xmz**2-s+kk)

c        IF(sn.ge.4d0.AND.z.gt.x)then   
c        If(s.ge.4*xmz**2.AND.t.ge.tlval.AND.t.le.thval) THEN
c                include './amplitudes/bornamp.h'
c                sigma(1) = 3*Bf*born
                sigma(1) = Bf*8*(t*t+u*u)/s/s
        !ONE LOOP CONTRIBUTION
        
c        If(ORD.eq.1) then 
        If(ORD.eq.1.AND.born.ge.0d0) then 
c                include './amplitudes/oneLampHG.h'
                !g01 = oneLoop/born 
                asp = as*4d0*pi/(1 +(bt0)*as*dlog(xmuR**2/s))
                zeta2 = pi2/6.0d0

c                oneL = (2*Cf*(5*s**2*(-1 + 5*zeta2) 
c     .              + 2*s*t*(-8 + 25*zeta2) + 2*t**2*(-8 + 25*zeta2) - 
c     .         (s**2 + 6*s*t + 6*t**2)*(0.0d0) - (s**2 + 2*s*t 
c     .         + 2*t**2)*(-6*zeta2)))/(Pi*s**2)

                oneL = (s*s + 2*s*t+2*t*t)*(8.0d0-7.0d0*zeta2)/pi/s/s
                oneL = -2.0d0*oneL*CF*asp*BF
                sigma(2) = oneL*2.0d0/sigma(1)
c                write(*,*)'Sigma 2 =',sigma(2)

c                sigma(2) = oneL*asp*2d0*4d0/3d0/2d0/pi/born
        Elseif(ORD.eq.2.and.born.gt.0d0) then
        !TWO LOOP CONTRIBUTION
c                include './amplitudes/oneLampHG.h'
                asp = as*4d0*pi/(1 +(bt0)*as*dlog(xmuR**2/s))
                sigma(2) = oneL*asp*2d0*cf/2d0/pi
                sigma(2) = sigma(2)/born

                sigma(3) = 0d0                     

c        If(M11switch.eq.1)then
c         call  CPP_WRAPPER(sn,un,m0m0,m0m1,m1m1,m0m2a,m0m2b,m0m2c,m0m2d)
c                twoL11 = m1m1
c                sigma(3) = twoL11*cf2/2d0
c                sigma(2) = 2*(as*4d0*pi)**2/4d0/pi2*sigma(3)/born
c                sigma(4) = 0d0
c        else

c                include './amplitudes/twoL02aamp.h'
                ! <M(0)|M(2)> contribution
                sigma(3) = sigma(3) + twoL02a*cf*xnf

c               include './amplitudes/twoL02camp.h'
               sigma(3) = sigma(3) + twoL02c*cf2

c               include './amplitudes/twoL02damp.h'
               sigma(3) = sigma(3) + twoL02d*CF/3d0

c              !< g02 = TwoLoop/born >
               asp = as*4d0*pi
               sigma(3) = 2*(asp)**2/4d0/pi2*sigma(3)/born
               sigma(2) = sigma(2) + sigma(3)
               
c               include './amplitudes/twoL02bamp.h'
               sigma(4) = twoL02b*CF
c               !< MULTIPLYING NFzfactor >
               sigma(4) = sigma(4)*(3d0*gvd**2+2d0*gvu**2+5d0*gau**2)
c               !< g02 >
               sigma(4) = 2*(asp)**2/4d0/pi2*sigma(4)/born
c       endif

        Else  ! if order is not 2 or 1 i.e. Leading order 
                sigma(2) = 0
                sigma(3) = 0
                sigma(4) = 0
        Endif

c        Else               ! if s and t are not physical
c                Sigma(1) = 0
c                Sigma(2) = 0
c                Sigma(3) = 0
c                Sigma(4) = 0
c        Endif
              return
        END

      function xlambda(ss,hh,zz)
      implicit double precision (a-h,o-z)
      ss2 = ss*ss
      zz2 = zz*zz
      hh2 = hh*hh
      xlambda = ss2+hh2+zz2-2d0*(ss*hh+ss*zz+zz*hh)
      return
      end
