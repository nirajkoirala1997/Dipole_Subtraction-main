c-------------------------------------------c
c-     2 -> 2 matrix element squared       -c
c-------------------------------------------c      

c     SM AND BSM
      subroutine setsigma(s,t,u,sigma)

        use handyG

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

        call clearcache

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
        Bf = 4d0/9d0*ALEM2*Pi**2
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
        If(s.ge.4*xmz**2.AND.t.ge.tlval.AND.t.le.thval) THEN
                include './amplitudes/bornamp.h'
                sigma(1) = 3*Bf*born
        !ONE LOOP CONTRIBUTION
        
c        If(ORD.eq.1) then 
        !TWO LOOP CONTRIBUTION
                sigma(3) = 0d0                     

        call  CPP_WRAPPER(sn,un,m0m0,m0m1,m1m1,m0m2a,m0m2b,m0m2c,m0m2d)
                twoL11 = m1m1
                sigma(3) = twoL11*cf2/2d0
                sigma(3) = 2*(as*4d0*pi)**2/4d0/pi2*sigma(3)/born
                sigma(2) = sigma(3)
                sigma(4) = 0d0
        Else               ! if s and t are not physical
                Sigma(1) = 0
                Sigma(2) = 0
                Sigma(3) = 0
                Sigma(4) = 0
        Endif
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
