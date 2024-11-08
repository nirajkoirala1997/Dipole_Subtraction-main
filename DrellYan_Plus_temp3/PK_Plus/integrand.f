      function flo2_Plus(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      double precision Kbar_qq_Plus,Ktil_qq_Plus
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq

      xa     = yy(1)
      xb     = yy(2)
      xc     = yy(3)

      delta  = 1d-6
c      xmin   = 0.0d0
c      xmax   = 1.0d0 - delta
c      xjac4   = (xmax-xmin)
c      x      = xmin+ xjac4*yy(4)
c
      almin = delta
      almax = 1.0d0
      al = almin*(almax/almin)**yy(4)
      aljac = al*dlog(almax/almin)
      xjac4 = aljac
      x = 1.0d0-al


      xjact = 2.0d0
      xjac = xjact*xjac4

      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps


        flo2_Plus  = 0.0d0
        PKplus_x = 0.0d0
        PKplus_1 = 0.0d0
        PKplus = 0.0d0
        PKReg    = 0.0d0
        PKDel    = 0.0d0

        sig1 = 0.0d0
        sig2 = 0.0d0
        sig3 = 0.0d0
        sig4 = 0.0d0

        do k = 1,2

c        do iplus = 0,1
c        if (x .le. 1d0-delta_cut) iplus = 0
c        if (x .ge. 1d0-delta_cut) iplus = 1

c        if (iplus .eq. 1) then

         if (k .eq. 1) call kinvar2_PK(x*xa,xb,xc,Qmass,p1,p2,p3,p4)
         if (k .eq. 2) call kinvar2_PK(xa,x*xb,xc,Qmass,p1,p2,p3,p4)

c        elseif (iplus .eq. 0) then
c        call kinvar2_PK(xa,xb,xc,Qmass,p1,p2,p3,p4)
c        endif

         scale = Qmass

        if (scale .ge. xlow .and. scale .le. xhigh)  then


         coef = Born_uU2eE(0,p1,p2,p3,p4)

         sp   = 2.0d0*dot(p1,p2)
         s12  = 2.0d0*dot(p1,p2)
         rsp  = dsqrt(sp) 
         pin  = 0.5d0*rsp
         pf   = 0.5d0*rsp
         flux = 4.0d0*pin*rsp
         
             xmuf = scale
             xmur = xmuf
            xmuf2 = xmuf*xmuf 
               AL = alphasPDF(xmur) 
              ALP = AL/2d0/Pi

            azmth = 2.0d0*pi
              ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

c            call getPKPlus(iplus,x,xmuf,p1,p2,p3,p4,SumPlus)
c           xone_plus_x2 = 1d0 + x**2
c           xone_minusx  = 1d0 - x
c           SumPlus = xone_plus_x2/xone_minus_x
c            SumPlus = (1d0+x**2)/(1d0-x) 

c        Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
              P_qq_Plus = Cf*(-1.0d0)*dlog(xmuf2/s12)*(1d0+x**2)/(1d0-x)
            aKbar_qq_Plus =Cf* 2d0*dLog(1d0 - x)/(1d0 - x) 
            aKtil_qq_Plus =Cf* 2d0*dLog(1d0-x)/(1d0 - x) 

c           SumPlus = P_qq_Plus + aKbar_qq_Plus + aKtil_qq_Plus
           SumPlus =  aKbar_qq_Plus 

            Cf = 4d0/3d0 
            sig = xl(1)*SumPlus*coef
c            if (iplus .eq. 0) sig3 = xl(1)*SumPlus*coef
  
c          if (iplus .eq. 1) then
            wgt = sig/flux*ps2*xjac*vwgt
            PKplus = xnorm*wgt/vwgt/2d0/eps

c          elseif (iplus .eq.0) then
c            wgt3 = sig3/flux*ps2*xjac*vwgt
c            PKplus_1= xnorm*wgt3/vwgt/2d0/eps

c         endif

           endif                ! bin choice
c           enddo

c         PK(k) = PKplus_x - PKplus_1 ! Plus is regulated at SumPlus
c Iplus will divide integral into two regions.
         PK(k) = PKplus
c	 print*,PKplus
         enddo

        flo2_Plus = ALP * ( PK(1) + PK(2) )
      return
      end
c---------------------------------------------------------------------
