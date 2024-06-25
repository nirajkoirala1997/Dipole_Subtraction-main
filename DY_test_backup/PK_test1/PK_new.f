! a -> ai+ i ||  q -> q g || P{a,ai} -> P{q,q} || a- Incoming
!                                               | ai- Undergoing Born
!                                               |    process               
! Kb = Kbar      
c----------------------------------------------------------- 
      subroutine getPKPlus(x,xmuf,p,xp1,xp2,SumP)
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

c        coef1 = (-1.0d0)*dlog(xmuf2/s12)      !{a,ai,i} => {q,q,g}
        Born = Born_uU2eE(0,p1,p2,p3,p4)
c        coef1 = coef1*Born1 
        coef = Born 

        ALLP(k) = Alp*(PqqP(x)+AKbarP_qq(x)+AKtilP_qq(x))*coef
c        ALLP(k) = Alp*PqqP(x)*coef


       enddo

       SumP(1) = AllP(1)+AllP(2)

      return
      end
cc----------------------------------------------------------- !Misc-functions
      subroutine getPKReg(x,xmuf,p,xp1,xp2,SumP)
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

c        coef1 = (-1.0d0)*dlog(xmuf2/s12)      !{a,ai,i} => {q,q,g}
        Born = Born_uU2eE(0,p1,p2,p3,p4)
c        coef1 = coef1*Born1 
        coef = Born 

       Areg = (AKbarreg_qq(x)+AKtilreg_qq(x))*coef
       AllP(k)= Alp*Areg

       enddo

       SumP(1) = AllP(1)+AllP(2)

      return
      end
cc----------------------------------------------------------- !Misc-functions
       subroutine getPKDel(x,xmuf,p,xp1,xp2,SumP)
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

c        coef1 = (-1.0d0)*dlog(xmuf2/s12)      !{a,ai,i} => {q,q,g}
        Born = Born_uU2eE(0,p1,p2,p3,p4)
c        coef1 = coef1*Born1 
        coef = Born 

       ADel = (AKbarD_qq(x)+AKtilD_qq(x))*coef
c       ADel = AKbarD_qq(x)*coef
c       ADel = AKtilD_qq(x)*coef
       AllP(k)= Alp*ADel
c       AllP(k)= Alp*coef
c       print*,AKbarD_qq(x)+AKtilD_qq(x)

       enddo

       SumP(1) = AllP(1)+AllP(2)

      return
      end
cc----------------------------------------------------------- !Misc-functions
                                                            
      double precision function PqqP(x)
      implicit double precision (a-h,o-z)
      data PI/3.141592653589793238462643D0/

      CF = 4.0d0/3.0d0
      PqqP = CF*(1.0d0+x*x)/(1.0d0-x)

c      XM=1.0D0-X
c      DLOGXM= DLOG(XM)
c      XD1=DLOGXM/XM
c      CQQB1PLUS=16.0D0*XD1
c      PqqP=CF*CQQB1PLUS
      return
      end

      double precision function Pqqreg(x)
      implicit double precision (a-h,o-z)
      CF = 4.0d0/3.0d0
      pqqreg = -CF*(1.0d0+x)
      return
      end

      double precision function AKbarP_qq(x)
      implicit double precision (a-h,o-z)
      
      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      AKbarP_qq = CF*2.0d0/(1.0d0-x)*dlgx 
      return
      end

      double precision function AKbarreg_qq(x)
      implicit double precision (a-h,o-z)
      
      dlgx = dlog((1.0d0-x)/x)
      CF = 4.0d0/3.0d0
      AKbarreg_qq = -CF*(1.0d0+x)*dlgx + (1-x) 
      return
      end

      double precision function AKbarD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      
      CF = 4.0d0/3.0d0
      AKbarD_qq = -CF*(5.0d0-pi*pi)
      return
      end

      double precision function AKtilP_qq(x)
      implicit double precision (a-h,o-z)
      
      CF = 4.0d0/3.0d0
      dlgx = dlog(1.0d0-x)
      AKtilP_qq = CF*2.0d0/(1.0d0-x)*dlgx 
      return
      end


      double precision function AKtilreg_qq(x)
      implicit double precision (a-h,o-z)
      external Pqqreg
      
      dlgx = dlog(1.0d0-x)
      AKtilreg_qq = Pqqreg(x)*dlgx
      return
      end

      double precision function AKtilD_qq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      CF = 4.0d0/3.0d0
      
      AKtilD_qq = -CF*pi*pi/3.0d0
      return
      end


      double precision function Pgq_reg(x) !same as Pgqb(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)

      Tr = 0.5d0
      Pgq_reg = Tr*(x**2 + (1-x)**2)

      return
      end

      double precision function aKbar_gq(x)
      implicit double precision (a-h,o-z)
      external Pgq_reg

      Tr = 0.5d0
      dlgx = dlog((1d0-x)/x)
      aKbar_gq = Pgq_reg(x)*dlgx +Tr*2d0*x*(1d0-x)
      return
      end


      double precision function aKtil_gq(x)
      implicit double precision (a-h,o-z)
      external Pgq_reg

      dlgx = dlog(1d0 - x)
      aKtil_gq = Pgq_reg(x)*dlgx 

      return
      end
