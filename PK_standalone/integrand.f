      function flo2_PK(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),zz(10)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),xp1(0:3),xp2(0:3)
     .          ,p(0:3,1:4),Born(1:2),all_PS(0:3,1:4,0:2)
      dimension SumP(1:2),SumK(1:2)
      dimension pa1(0:3),pa2(0:3),pa3(0:3),pa4(0:3)
      dimension pb1(0:3),pb2(0:3),pb3(0:3),pb4(0:3)
      dimension xinvmass(0:2)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      character*50 name
      common/energy/s
      common/pdfname/name
      common/factscale/xmuf
      common/leg_choice/leg
      common/usedalpha/AL,ge
      common/distribution/xq
      external PqqP,Pqqreg,Born_uU2eE
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg

      xa     = yy(1)
      xb     = yy(2)
      xc     = yy(c)
      xmin   = 0.0d0
      xmax   = 1.0d0 - 1d-3
      xjac   = (xmax-xmin)
      x      = xmin+ xjac*yy(4)
      yy(4)  = x

      sp     = xa*xb*s
      rsp    = dsqrt(sp) 
      xjac = xjac*2.0d0
      xnorm=hbarc2

      eps = 1.0d0
      xlow = xq - eps
      xhigh = xq + eps


c     all_PS -> (4momenta,particleID,kinematic type)
c        call kinvar2_PK(yy,xinvmass,all_PS)

        flo2_PK  = 0.0d0
        PKplus_x = 0.0d0
        PKplus_1 = 0.0d0
        PKReg    = 0.0d0
        PKDel    = 0.0d0

        sig1 = 0.0d0
        sig2 = 0.0d0
        sig3 = 0.0d0
        sig4 = 0.0d0

        do k = 1,2

        do iplus = 0,1

        if (iplus .eq. 1) then

        if (k .eq. 1) then
        call kinvar2_PK(x*xa,xb,xc,xxinvmass,p1,p2,p3,p4)
        elseif (k .eq. 2) then
        call kinvar2_PK(xa,x*xb,xc,xxinvmass,p1,p2,p3,p4)
        endif

        elseif (iplus .eq. 0) then
        call kinvar2_PK(xa,xb,xc,xxinvmass,p1,p2,p3,p4)
        endif

        scale = xxinvmass

        pin  = 0.5d0*rsp
        flux_1 = 4.0d0*pin*rsp
        flux_x = 4.0d0*x*pin*rsp
        coef = Born_uU2eE(0,p1,p2,p3,p4)


        if (scale .ge. xlow .and. scale .le. xhigh) then   

             xmuf = scale
             xmur = xmuf
            xmuf2 = xmuf*xmuf 
c               AL = alphasPDF(xmur)
c              ALP = AL/2d0/Pi
              s12 = 2d0*dot(p1,p2)

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)


            call getPK(x,xmuf,p1,p2,p3,p4,SPKplus,SPKReg,SPK1,SPKDel)
                
c           coefPterm   = (-1.0d0)*dlog(xmuf2/s12/x)
c           SumP(1) = coefPterm*PqqP(x) 

c           SumK(1) = AKbarP_qq(x)+AKtilP_qq(x)

c            sig1 = xl(1)*( SumP(1) + SumK(1) ) 

            if (iplus .eq. 1) then
            sig1 = xl(1)*(SumP + SumK) 
            sig2 = xl(1)*Reg
            elseif (iplus .eq. 0) then
            sig3 = xl(1)*Plus_1
            sig4 = xl(1)*Del
            endif
  
            sig = Alp*sig1*coef

            azmth = 2.0d0*pi
            pf = 0.5d0*rsp
            ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth


            if (iplus .eq. 1) then
            wgt1 = sig1/flux_x*ps2*xjac*vwgt
            wgt2 = sig2/flux_x*ps2*xjac*vwgt
            PKplus_x(k) = xnorm*wgt1/vwgt/2d0/eps
            PKReg(k) = xnorm*wgt2/vwgt/2d0/eps

           elseif (iplus .eq.0) then
            wgt3 = sig3/flux_x*ps2*xjac*vwgt
            wgt4 = sig4/flux_x*ps2*xjac*vwgt
            PKplus_1(k)= xnorm*wgt3/vwgt/2d0/eps
            PKDel(k) = xnorm*wgt4/vwgt/2d0/eps
            endif
            enddo

            PK(k) = PKplus_x - PKplus1 + PKReg + PKDel
            enddo

            flo2_PK = PK1(1) + PK(2)

      return
      end

c---------------------------------------------------------------------
      subroutine momenta_giver(kin_type,all_PS,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
c     all_PS(4momenta,particleID,kinematic type)
       dimension all_PS(0:3,1:4,0:2)
       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
         do j=0,3
           p1(j) = all_PS(j,1,kin_type) 
           p2(j) = all_PS(j,2,kin_type) 
           p3(j) = all_PS(j,3,kin_type) 
           p4(j) = all_PS(j,4,kin_type) 
         enddo
      end
c---------------------------------------------------------------------
 

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
