      function flo2_PK(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2)
      dimension SumP(1:2),SumK(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/leg_choice/leg
      common/usedalpha/AL,ge
      common/distribution/xq

      xa     = yy(1)
      xb     = yy(2)
      delta  = 1d-5

      xmin   = 0.0d0
      xmax   = 1.0d0 - delta
      xjac   = (xmax-xmin)
      x      = xmin+ xjac*yy(4)

c      almin = delta
c      almax = 1.0d0
c      al = almin*(almax/almin)**yy(4)
c      aljac = al*dlog(almax/almin)
c      aljac = -dlog(delta)
c      xjac = aljac
c      x = 1.0d0-al

      sp     = xa*xb*s
      rsp    = dsqrt(sp) 
      xjac = xjac
      xnorm=hbarc2

      eps = 1.0d0
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
c        write(*,*)'Flux =', flux_x

        scale2 = 2.0d0*dot(p1,p2)
        scale = dsqrt(scale2)

      flo2_PK  = 0.0d0
      PKplus_x = 0.0d0
      PKplus_1 = 0.0d0
      PKRegDel = 0.0d0
      PKplus   = 0.0d0

c        goto 101

        if (scalex .ge. xlow .and. scalex .le. xhigh)  then   

            xmuf = scalex
            xmur = scalex
c            AL = alphasPDF(xmur)
            AS = 1.0d0
            ALP = AS

            call getPK(2,x,xmuf,p,xp1,xp2,SumP)
                
            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(1)* SumP(1)  !  [qq lum]

c        write(*,*)'Flux =', scalex
            sig = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt_x = sig/flux_x*ps2*xjac*vwgt
           PKplus_x = xnorm*wgt_x/vwgt/2d0/eps

        write(*,*)'PKplus_x =', SumP(1) 
         endif
c        write(*,*)'PKplus_x out =', PKplus_x 

c 101    continue

c        goto 102

        if (scale .ge. xlow .and. scale .le. xhigh)  then   

            xmuf = scale
            xmur = scale
c            AL = alphasPDF(xmur)
            AS = 1.0d0
            ALP = AS

            call getPK(2,x,xmuf,p,xp1,xp2,SumP)
                
            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = xl(1)* SumP(1)  !  [qq lum]

            sig = Alp*sig1

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt = sig/flux_1*ps2*xjac*vwgt
            PKplus_1 = xnorm*wgt/vwgt/2d0/eps

         endif

c  102     continue


c         PKplus = PKplus_x - PKplus_1

c        write(*,*)'PKplus =', PKplus_x, PKplus_1

c        if (scalex .ge. xlow .and. scalex .le. xhigh) then
c
c            xmuf = scalex
c            xmur = scalex
cc            AL = alphasPDF(xmur)
c            Alp = 1.0d0
c
c            call getPK(1,x,xmuf,p,xp1,xp2,SumP)
cc                
c            call pdf(xa,xmuf,f1)
c            call pdf(xb,xmuf,f2)
c            call setlum(f1,f2,xl)
c
c            sig1 = xl(1)* (SumP(1))  !  [qq lum]
c
c            sig = Alp*sig1
c
c            azmth = 2.0d0*pi
c            pf = 0.5d0*rsp
c            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth
c
c            wgt = sig/flux_x*ps2*xjac*vwgt
c            PKreg = xnorm*wgt/vwgt/2d0/eps
c
c         endif
c
c
c         if (scale .ge. xlow .and. scale .le. xhigh)  then   
c
c            xmuf = scale
c            xmur = scale
c            AL = alphasPDF(xmur)
c            AS = 1.0d0
c            ALP = AS
c
c            call getPK(0,x,xmuf,p,xp1,xp2,SumP)
c                
c            call pdf(xa,xmuf,f1)
c            call pdf(xb,xmuf,f2)
c            call setlum(f1,f2,xl)
c
cc            sig1 = xl(1)* (SumP(1)+SumK(1))  !  [qq lum]
c            sig1 = xl(1)* SumP(1)  !  [qq lum]
c
c            sig = Alp*sig1
c
c            azmth = 2.0d0*pi
c            pf = 0.5d0*rsp
c            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth
c
c            wgt = sig/flux_1*ps2*xjac*vwgt
c            PKdel = xnorm*wgt/vwgt/2d0/eps
c
c         endif

c         flo2_PK = PKplus + PKreg + PKdel
         flo2_PK = PKplus_x

     return
      end

c---------------------------------------------------------------------
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
