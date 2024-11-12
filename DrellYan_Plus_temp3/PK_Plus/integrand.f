      function flo2_Plus(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
      double precision Kbar_qq_Plus,Ktil_qq_Plus
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/plus_cutoff/delta

      xa     = yy(1)
      xb     = yy(2)
      xc     = yy(3)
c
c      delta  = 1d-6
c      xmin   = 0.0d0
c      xmax   = 1.0d0 - delta
c      xjac4   = (xmax-xmin)
c      x      = xmin+ xjac4*yy(4)
c
c      delta = 1d-6
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


         PKplus = 0.0d0
        SumPlus = 0.0d0

        do k = 1,2

         if (k .eq. 1) call kinvar2_PK(x*xa,xb,xc,Qmass,p1,p2,p3,p4)
         if (k .eq. 2) call kinvar2_PK(xa,x*xb,xc,Qmass,p1,p2,p3,p4)

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

                P_qq_Plus = PqqP(x)
                P_qq_Plus = P_qq_Plus*(-1.0d0)*dlog(xmuf2/s12)
            aKbar_qq_Plus = aKbarP_qq(x) 
            aKtil_qq_Plus = aKtilP_qq(x) 

c           SumPlus = P_qq_Plus + aKbar_qq_Plus + aKtil_qq_Plus
c           SumPlus = aKtil_qq_Plus
           SumPlus = aKbar_qq_Plus
c           SumPlus = P_qq_Plus

            sig = xl(1)*SumPlus*coef
            wgt = sig/flux*ps2*xjac*vwgt
            PKplus = xnorm*wgt/vwgt/2d0/eps

           endif                ! bin choice
          PK(k) = PKplus
         enddo

        flo2_Plus = ALP * ( PK(1) + PK(2) )
      return
      end

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      function flo2_PKDel(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     .         ,f1(-6:6),f2(-6:6),xl(15)
     .         ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     .         ,p(0:3,1:4),Born(1:2)
     .         ,PK(1:2)
      dimension AllPK(1:2),AllK(1:4)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/energy/s
      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/plus_cutoff/delta

c      xa     = yy(1)
c      xb     = yy(2)

c      delta = 1d-3
      TAUH = xq*xq/s

      xa = (1D0-TAUH)*yy(1)    + TAUH
c      xb = (1D0-TAUH/xa)*yy(2) + TAUH/xa
      xb = TAUH/xa
      xc = yy(2)

      sp     = xa*xb*s
      rsp    = dsqrt(sp)
      xjac = 2.0d0*(1.0d0-TAUH)
      xnorm=hbarc2*2.0d0*xq/s/xa

        call kinvar2_PK(xa,xb,xc,xinvmass,p1,p2,p3,p4)

        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp
        s12 = 2d0*dot(p1,p2)

        scale = xq

      flo2_PKDel  = 0.0d0
      sig1 = 0d0
      Pqq_int_ = 0d0


            xmuf = scale
            xmuf2 = xmuf*xmuf
            xmur = scale
            AL = alphasPDF(xmur)
            ALP = AL/2d0/Pi
             Cf = 4d0/3d0

         do k=1,2  ! for two Legs I can multiply by 2 also instead of using loop.
         ! corrspond to g(1)*(theta(1-delta) -theta(0) 


c       ~~~~[ P term ]~~~~c
c          Pqq_int_= Pqq_int_P(1d0-delta) - Pqq_int_P(0d0)
          Pqq_int_= -3d0/2d0-0.5d0*(delta-4d0)*delta - 2d0*dLog(delta) 
          Pqq_int_= (-1.0d0)*dlog(xmuf2/s12)*Pqq_int_*Cf 

c the coefficients log(muf/s12) is there as overall in eq.(10.25)

c       ~~~~[ Kb term ]~~~~c
          aKb_qq_int_ = Cf*Pi**2/3d0 - Cf*dLog(delta)**2 ! -2PolyLog[2,delta] !! Here Ignoring PolyLog as delta is very small. 


c       ~~~~[ Kt term ]~~~~c
          aKt_qq_int_ = -Cf*dLog(delta)**2

c            PK(k) = Pqq_int_ + aKb_qq_int_ + aKt_qq_int_ 
c            PK(k) = aKt_qq_int_ 
          PK(k) = aKb_qq_int_
c          PK(k) = Pqq_int_*Cf 
         enddo
        coef = Born_uU2eE(0,p1,p2,p3,p4)

         SumDel = PK(1)+PK(2)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)

            sig1 = ALP * xl(1) * coef * SumDel  

c          [Note-] 
c          qq lum and overall factor al/2/pi and Born with x=1
c          kinematics which is component of g(1) in our notation.


            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

            wgt1 = sig1/flux_1*ps2*xjac*vwgt
            PKDel = xnorm*wgt1/vwgt

         flo2_PKDel = PKDel

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms .end]
