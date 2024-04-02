c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Plus Terms]
      function flo2_PK(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq

      xa     = yy(1)
      xb     = yy(2)
      delta  = 1d-5

c      xmin   = 0.0d0
c      xmax   = 1.0d0 - delta
c      xjac   = (xmax-xmin)
c      x      = xmin+ xjac*yy(4)

      almin = delta
      almax = 1.0d0
      al = almin*(almax/almin)**yy(4)
      aljac = al*dlog(almax/almin)
      xjac = aljac
      x = 1.0d0-al

      sp     = xa*xb*s
      rsp    = dsqrt(sp) 
      xjac = xjac
      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)

            do i=0,3
             xp1(i) = x*p1(i)
             xp2(i) = x*p2(i)
            enddo

        call p1d_to_p2d_4(p1,p2,p3,p4,p)

        scalex2 = 2.0d0*x*dot(p1,p2)
        scalex = dsqrt(scalex2)
        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp
        flux_x = 4.0d0*x*pin*rsp

        scale2 = 2.0d0*dot(p1,p2)
        scale = dsqrt(scale2)

      flo2_PK  = 0.0d0
      PKplus_x = 0.0d0
      PKplus_1 = 0.0d0
      PKplus   = 0.0d0
      PKterms  = 0.0d0
      sig1 = 0d0

        if (scalex .ge. xlow .and. scalex .le. xhigh)  then   

            xmuf = scalex
            xmur = scalex
            AL = alphasPDF(xmur)
            AS = 1.0d0
            ALP = AS

            call getPKPlus(x,xmuf,p,xp1,xp2,SumPlus)
                
            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

c            sig1 = xl(1)* (SumP(1)+SumK(1))  !  [qq lum]
            sig1 = xl(1)* SumPlus  !  [qq lum]
            sig = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt_x = sig/flux_x*ps2*xjac*vwgt
           PKplus_x = xnorm*wgt_x/vwgt/2d0/eps

         endif

        if (scale .ge. xlow .and. scale .le. xhigh)  then   

            xmuf = scale
            xmur = scale
c           AL = alphasPDF(xmur)
            AS = 1.0d0
            ALP = AS

            call getPKPlus(x,xmuf,p,xp1,xp2,SumPlus)
                
            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(1)* SumPlus  !  [qq lum]
            sig1 = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt1 = sig1/flux_1*ps2*xjac*vwgt
            PKplus_1 = xnorm*wgt1/vwgt/2d0/eps

         endif

         PKterms = PKplus_x - PKplus_1
         flo2_PK = PKterms

      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      function flo2_PKReg(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2)
      dimension SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq

      xa     = yy(1)
      xb     = yy(2)
      x      = yy(4)

      sp     = xa*xb*s
      rsp    = dsqrt(sp) 
      xjac   = 1d0
      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)

            do i=0,3
             xp1(i) = x*p1(i)
             xp2(i) = x*p2(i)
            enddo

        call p1d_to_p2d_4(p1,p2,p3,p4,p)

        scalex2 = 2.0d0*x*dot(p1,p2)
        scalex = dsqrt(scalex2)
        pin  = 0.5d0*rsp
        flux_x = 4.0d0*x*pin*rsp

        PKReg = 0.0d0


        if (scalex .ge. xlow .and. scalex .le. xhigh) then   

            xmuf = scalex
            xmur = scalex
            AL = alphasPDF(xmur)
            AS = 1.0d0
            ALP = AS

            call getPKReg(x,xmuf,p,xp1,xp2,SumReg)
                
            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(1)* SumReg  !  [qq lum]

            sig = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt_x = sig/flux_x*ps2*xjac*vwgt
           PKReg = xnorm*wgt_x/vwgt/2d0/eps

         endif
         flo2_PKReg = PKReg

      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      function flo2_PKDel(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     .         ,f1(-6:6),f2(-6:6),xl(15)
     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .         ,p(0:3,1:4),Born(1:2)
     .         ,SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq

      xa     = yy(1)
      xb     = yy(2)

      sp     = xa*xb*s
      rsp    = dsqrt(sp) 
      xjac = 1.0d0
      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps

        call kinvar2(yy,xinvmass,p1,p2,p3,p4)

        call p1d_to_p2d_4(p1,p2,p3,p4,p)

        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp

        scale2 = 2.0d0*dot(p1,p2)
        scale = dsqrt(scale2)

      flo2_PKDel  = 0.0d0
      sig1 = 0d0


        if (scale .ge. xlow .and. scale .le. xhigh)  then   

            xmuf = scale
            xmur = scale
            AL = alphasPDF(xmur)
            AS = 1.0d0
            ALP = AS

            call getPKDel(1.0d0,xmuf,p,p1,p2,SumDel)
                
            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(1)* SumDel  !  [qq lum]

            sig1 = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt1 = sig1/flux_1*ps2*xjac*vwgt
            PKDel = xnorm*wgt1/vwgt/2d0/eps

         endif
         flo2_PKDel = PKDel

      return
      end


c--------------------------------------------------------------------o
c     [u U -> e E]  Born 
c--------------------------------------------------------------------o
       function Born_uU2eE(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       common/usedalpha/AL,ge
c       ge=1d0/128d0
       e= DSQRT(ge*4.d0*PI)
      IF(k .eq. 0)  CF =  1d0                   !Leading Order K=0
      IF(k .eq. 1)  CF = -4d0/3d0               !leg 1
      IF(k .eq. 2)  CF = -4d0/3d0               !Leg 2
      s13 =  2.0d0*dot(p1,p3) ! t
      s23 =  2.0d0*dot(p1,p4) ! u
      s12 =  2.0d0*dot(p1,p2) ! s
      qu2 = 1d0!4d0/9d0

      Born_uU2eE= CF*(2*e**4*qu2*(-2*s13*s23 + s12*(s13 +
     .            s23)))/(3d0*s12**2)
       return
       end
c---------------------------------------------------------------------
