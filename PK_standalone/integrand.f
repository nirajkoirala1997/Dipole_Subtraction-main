      function flo2_PK(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6),xl(1)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/leg_choice/leg
      common/usedalpha/AL
      common/distribution/xq
      


      xa     = yy(1)
      xb     = yy(2)
       x     = yy(4)
      rsp = dsqrt(xa*xb*s)

      ipass1 = 0

        eps = 0.5d0*2d0
       xlow = xq - eps
      xhigh = xq + eps

      xcut = xq - 10.0d0

      if (rsp .gt. xcut) ipass1 = 1
      if (ipass1 .eq.1) then

      call kinvar2(yy,xinvmass,p1,p2,p3,p4)
      call p1d_to_p2d_4(p1,p2,p3,p4,p)

      scale  = xinvmass
c      if( scale .gt. 100d0) then
        ipass = 0
        flo2_PK = 0
        if ( scale .ge. xlow .and. scale .le. xhigh) ipass=1
         if ( ipass .eq. 1 ) then
            xmuf=scale
            xmur=scale

       call pdf(xa,xmuf,f1)
       call pdf(xb,xmuf,f2)
       call setlum(f1,f2,xl)
       AL = alphasPDF(xmur)
      if (leg .eq. 1) call PKterm1(p,x,SumP,SumK)
      if (leg .eq. 2) call PKterm2(p,x,SumP,SumK)
c      call getPK(leg,x,xmuf,p,SumP,SumK)

        do i=0,3
        xp1(i) = x*p1(i)
        xp2(i) = x*p2(i)
        enddo

c      if (leg .eq. 1) Born(leg)=Born_uU2eE(0,xp1,p2,p3,p4)
c      if (leg .eq. 2) Born(leg)=Born_uU2eE(0,p1,xp2,p3,p4)
      if (leg .eq. 1) Born(leg)=Born_uU2eE(0,xp1,p2,p3,p4)
      if (leg .eq. 2) Born(leg)=Born_uU2eE(0,p1,xp2,p3,p4)

       PK =SumP+SumK

       sig = xl(1)*Born(leg)* PK

c            xnorm=hbarc2/16d0/pi/(s*xa*xb)
       xnorm=hbarc2/16d0/pi/s
       wgt=xnorm*sig*vwgt
       flo2_PK=wgt/vwgt
       return              ! There is a factor of 1/2
         else                   ! from xeps interval
            flo2_PK=0d0
            return
         endif
        endif

         return
      end

c---------------------------------------------------------------------
c     [u U -> e E]  Born 
c--------------------------------------------------------------------o
       function Born_uU2eE(k,p1,p2,p3,p4)
       implicit double precision (a-h,o-z)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       parameter(PI=3.141592653589793238D0)
       common/usedalpha/AL
       ge=0.007547169811320755d0
c       Al=0.118d0
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
c---------------------------------------------------------------------

