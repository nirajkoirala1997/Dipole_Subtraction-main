cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Plus Terms]
c      subroutine getPKPlus(iplus,x,xmuf,p1,p2,p3,p4,SumPlus)
c      implicit double precision (a-h,o-z)
c      parameter (pi=3.14159265358979d0)
c      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3)
c      common /usedalpha/ AL,ge 
c      external PqqP,Pqqreg
c      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
c      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
c      external aKbar_gq,aKtil_gq,Pgq_reg 
c
c        s12 = 2d0*dot(p1,p2)
c        Cf = 4d0/3d0                      
c        Tr = 0.5d0
c        Alp = Al/2.0d0/pi
c        Alp = 1.0d0
c        xmuf2 = xmuf*xmuf
c
c          Pplus = 0.0d0
c       SumPlus  = 0.0d0
c         SumReg = 0.0d0
c         SumDel = 0.0d0
c
c        if ( iplus .eq. 1 ) then
c          Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)
cc          write(*,*)'iplus, s12 =', iplus, s12
c        SumPlus = Pplus + AKbarP_qq(x) + AKtilP_qq(x)
c
c        elseif( iplus .eq. 0 ) then
cc          write(*,*)'iplus, s12 =', iplus, s12
c          Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
c        SumPlus = Pplus + AKbarP_qq(x) + AKtilP_qq(x)
c        endif
c
c      return
c      end
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
c      subroutine getPKReg(x,xmuf,p1,p2,p3,p4,SumReg)
c      implicit double precision (a-h,o-z)
c      parameter (pi=3.14159265358979d0)
c      dimension AllP(1:4),AllK(1:4),SumP(1:2),SumK(1:2)
c      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
c      dimension  xp1(0:3),xp2(0:3)
c      common /usedalpha/ AL,ge 
c      external Born_uU2eE
c      external PqqP,Pqqreg
c      external AKbarP_qq,AKbarreg_qq,AKbarD_qq
c      external AKtilP_qq,AKtilreg_qq,AKtilD_qq
c      external aKbar_gq,aKtil_gq,Pgq_reg
c
c       SumP(1) = 0d0
c       SumK(1) = 0d0
c       SumP(2) = 0d0
c       SumK(2) = 0d0
c
c        s12 = 2d0*dot(p1,p2)
c        Cf = 4d0/3d0
c        Tr = 0.5d0
cc        Alp = Al/2.0d0/pi
c        Alp = 1.0d0
c        xmuf2 = xmuf*xmuf
c
c        Born = Born_uU2eE(0,p1,p2,p3,p4)
c        coef = Born
c
c       Areg = Alp*(AKbarreg_qq(x)+AKtilreg_qq(x))*coef
c
c       SumReg = Areg
c      return
c      end
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
        delta = 1d-6
c        Alp = Al/2.0d0/pi
        xmuf2 = xmuf*xmuf

      do k = 1,2

        Born = Born_uU2eE(0,p1,p2,p3,p4)
        coef = Born

       P_qq_delta  = 2d0*dLog(delta) + 3d0/2d0  
       P_qq_delta = P_qq_delta * (-1.0d0)*dlog(xmuf2/s12)
       aKb_qq_delta = dLog(delta)**2 - 4d0*Pi**2/3d0 
       aKt_qq_delta = dLog(delta)**2 
c       del_all = Cf*(P_qq_delta + aKb_qq_delta + aKt_qq_delta) ! Cf for the functions are already multiplied
       del_all = Cf*aKb_qq_delta  ! Cf for the functions are already multiplied

c       AllP(k)= (AKbarD_qq(x)+AKtilD_qq(x) + del_all)*coef
       AllP(k)= (AKbarD_qq(x)+ del_all)*coef


c       AllP(k)= 2d0*dLog(1d-6) + 3d0/2d0  
c       AllP(k)= AllP(k)*Coef
      enddo

       SumDel = AllP(1)+AllP(2)

      return
      end
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Misc functions]
