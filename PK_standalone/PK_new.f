c---------------------------------------------------------------------------- 
      subroutine getPK(iplus,x,xmuf,p1,p2,p3,p4,SumPlus,SumReg,SumDel)
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
        Tr = 0.5d0
        Alp = Al/2.0d0/pi
        xmuf2 = xmuf*xmuf

          Pplus = 0.0d0
       SumPlus  = 0.0d0
         SumReg = 0.0d0
         SumDel = 0.0d0

        if ( iplus .eq. 1 ) then
          Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12/x)
        SumPlus = Pplus + AKbarP_qq(x) + AKtilP_qq(x)
         SumReg = (AKbarreg_qq(x)+AKtilreg_qq(x))
         SumPlus = PqqP(x)

        elseif( iplus .eq. 0 ) then
          Pplus = PqqP(x)*(-1.0d0)*dlog(xmuf2/s12)    ! here x=1
         SumDel = (AKbarD_qq(x)+AKtilD_qq(x))
        SumPlus = Pplus + AKbarP_qq(x) + AKtilP_qq(x)

        endif

      return
      end

c----------------------------------------------------------- !Misc-functions
                                                            
      double precision function PqqP(x)
      implicit double precision (a-h,o-z)
      data PI/3.141592653589793238462643D0/

      CF = 4.0d0/3.0d0
c      PqqP = CF*(1.0d0+x*x)/(1.0d0-x)

      XM=1.0D0-X
      LOGXM= DLOG(XM)
      XD1=LOGXM/XM
      CQQB1PLUS=16.0D0*XD1
      PqqP=CF*CQQB1PLUS
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
