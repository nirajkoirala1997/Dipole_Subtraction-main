c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Plus Terms]
      subroutine getPKPlus(iplus,x,xmuf,p1,p2,p3,p4,SumPlus)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common /usedalpha/ AL,ge 
      external PqqP,Pqqreg
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg 

        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0                      
        C_F = 4d0/3d0                      
        Tr = 0.5d0
        Alp = Al/2.0d0/pi
        xmuf2 = xmuf*xmuf
         pi2_3= Pi**2/3d0
         pi2  = Pi**2

          Pplus = 0.0d0
       SumPlus  = 0.0d0
         SumReg = 0.0d0
         SumDel = 0.0d0
c
c        if ( iplus .eq. 1 ) then
c
c        Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)
c        SumPlus = Pplus + AKtilP_qq(x) + AKbarP_qq(x)
cc        SumPlus = AKbarP_qq(x)
c
c        elseif( iplus .eq. 0 ) then
c
c        Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
c        SumPlus = Pplus + AKbarP_qq(x) + AKtilP_qq(x)
cc        SumPlus = AKbarP_qq(x)
c        endif
c~~~~~~~~~~~[ Delta functions Modified ]~~~~~~~~~~~
       if (iplus .eq. 1) then
         P_qq_delta = C_F * (3.0d0 / 2.0d0)
      else
         P_qq_delta = 0.0d0
      endif

c      P_gg_delta = (11.0d0 / 6.0d0 * C_A - 2.0d0 / 3.0d0 * N_f * T_R)

      if (iplus .eq. 0) then
         ! original definition of +-distribution
         Kbar_qq_delta = C_F * (-(5.0d0 - pi2))
      elseif (iplus .eq. 1) then
         ! alternative definition of +-distribution
         Kbar_qq_delta = C_F * (-(5.0d0 - pi2) - pi2_3)
      endif
         Kt_qq_delta = C_F * (-pi2 / 3.)
      if (iplus .eq. 1) then
         P_qq = -C_F * (1.0d0 + x)
      else
         P_qq = 0.0d0
      endif
c~~~~~~~~~~~[ Plus functions Modification ]~~~~~~~~~~~
c   P-terms Only
      iswitch_off_intXplus = 0
      if (iplus .eq. 1) then
         P_qq_plus = C_F * 2.0d0 / (1.0d0 - x)
      else
         P_qq_plus = C_F * (1.0d0 + x**2) / (1.0d0 - x)
      endif

      if (iswitch_off_intXplus .eq. 1) then
       xintP_qq_plus = 0.0d0
      else
         if (iplus .eq. 1) then
       xintP_qq_plus = -C_F * 2.0d0 * log(1.0d0 - x)
         else
       xintP_qq_plus = C_F*(-0.5d0 * x**2 - x - 2.0d0 * log(1.0d0 - x))
         endif
      endif
      SumPlus = P_qq_Plus + XintP_qq_plus 
c      print*,SumPlus,P_qq_Plus,XintP_qq_plus 
c   Kbar-terms Only

      if (iplus .eq. 0) then
         ! Original definition of +-distribution
         Kbar_qq = C_F * (-(1.0d0 + x)*log((1.0d0 - x) / x)+(1.0d0 - x))
      elseif (iplus .eq. 1) then
         Kbar_qq = C_F * (-(1.0d0 + x) * log((1.0d0 - x) / x) 
     .            + (1.0d0 - x) - 2.0d0 * log(x) / (1.0d0 - x))
      endif

      if (iplus .eq. 0) then
         Kbar_qq_plus = C_F*((2.0d0/(1.0d0 - x))*log((1.0d0 - x) / x))
      elseif (iplus .eq. 1) then
         Kbar_qq_plus = C_F * ((2.0d0 / (1.0d0 - x)) * log(1.0d0 - x))
      endif

      if (iswitch_off_intXplus .eq. 1) then
         xintKbar_qq_plus = 0.0d0
      else
         if (iplus .eq. 0) then
          xintKbar_qq_plus = C_F * (-(log(1.0d0 - x)**2)) 
c     .                      - 2.0d0 * gsl_sf_dilog(1.0d0 - x) + pi2_3)
         elseif (iplus .eq. 1) then
            xintKbar_qq_plus = C_F * (-(log(1.0d0 - x)**2))
         endif
      endif

c K_tilde-terms

c      Kt_qq = P_qq_reg(x) * log(1.0d0 - x)

      Kt_qq_plus = C_F * (2.0d0 / (1.0d0 - x)) * log(1.0d0 - x)

      if (iswitch_off_intXplus .eq. 1) then
         xintKt_qq_plus = 0.0d0
      else
         xintKt_qq_plus = -C_F * (log(1.0d0 - x)**2)
      endif





      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      subroutine getPKReg(x,xmuf,p1,p2,p3,p4,SumReg)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      dimension  xp1(0:3),xp2(0:3)
      common /usedalpha/ AL,ge 
      external Born_uU2eE
      external PqqP,Pqqreg
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg

       SumP(1) = 0d0
       SumK(1) = 0d0
       SumP(2) = 0d0
       SumK(2) = 0d0

        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0
        Tr = 0.5d0
c        Alp = Al/2.0d0/pi
        Alp = 1.0d0
        xmuf2 = xmuf*xmuf

        Born = Born_uU2eE(0,p1,p2,p3,p4)
        coef = Born

       Areg = Alp*(AKbarreg_qq(x)+AKtilreg_qq(x))*coef

       SumReg = Areg
      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Delta terms]
      subroutine getPKDel(x,xmuf,p,xp1,xp2,SumDel)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      dimension  xp1(0:3),xp2(0:3)
c      common /usedalpha/ AL,ge 
      external Born_uU2eE
      external PqqP,Pqqreg
      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
      external aKbar_gq,aKtil_gq,Pgq_reg

       SumP(1) = 0d0
       SumK(1) = 0d0
       SumP(2) = 0d0
       SumK(2) = 0d0

      call p2d_to_p1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0
        Tr = 0.5d0
c        Alp = Al/2.0d0/pi
        Alp = 1.0d0
        xmuf2 = xmuf*xmuf

      do k = 1,2

        Born = Born_uU2eE(0,p1,p2,p3,p4)
        coef = Born

       AllP(k)= (AKbarD_qq(x)+AKtilD_qq(x))*coef
      enddo

       SumDel = AllP(1)+AllP(2)

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ **end** ]
