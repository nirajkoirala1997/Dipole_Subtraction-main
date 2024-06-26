      function flo2_Plus(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10),PK(1:2)
      dimension f1(-6:6),f2(-6:6),xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension xp1(0:3),xp2(0:3)
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
      xmin   = 0.0d0
      xmax   = 1.0d0 - 1d-5
      xjac   = (xmax-xmin)
      x      = xmin+ xjac*yy(4)

c      sp     = xa*xb*s
c      rsp    = dsqrt(sp) 
      xnorm=hbarc2

      eps = 0.5d0
      xlow = xq - eps
      xhigh = xq + eps


        flo2_PK  = 0.0d0
        PKplus_x = 0.0d0
        PKplus_1 = 0.0d0
        PKReg    = 0.0d0
        PKDel    = 0.0d0

        sig1 = 0.0d0
        sig2 = 0.0d0
        sig3 = 0.0d0
        sig4 = 0.0d0

c        call kinvar2_PK(xa,xb,xc,p1,p2,p3,p4)

c        do i = 0,3
c        xp1(i) = x*p1(i)
c        xp2(i) = x*p2(i)
c        enddo


        do k = 1,2

        do iplus = 0,1

        if (iplus .eq. 1) then

         if (k .eq. 1) then
        call kinvar2_PK(x*xa,xb,xc,Qmass,p1,p2,p3,p4)
         elseif (k .eq. 2) then
         call kinvar2_PK(xa,x*xb,xc,Qmass,p1,p2,p3,p4)
         endif


        elseif (iplus .eq. 0) then
        call kinvar2_PK(xa,xb,xc,Qmass,p1,p2,p3,p4)
        endif

         scale = Qmass

        if (scale .ge. xlow .and. scale .le. xhigh) then   


c        if ( x .le. 0.5d0) then
c        if (iplus .eq. 0) then
c        write(*,*)'scale =', iplus, x,xa,scale
c        elseif (iplus .eq. 1) then
c        write(*,*)'scale =', iplus, x, xa,scale
c        write(*,*)' '
c        endif
c        endif

         coef = Born_uU2eE(0,p1,p2,p3,p4)

         sp   = 2.0d0*dot(p1,p2)
c        write(*,*)'xa =', iplus, xa,x
c        write(*,*) 'Ecm =', iplus, rsp
        
         rsp  = dsqrt(sp) 
         pin  = 0.5d0*rsp
         pf   = 0.5d0*rsp
         flux = 4.0d0*pin*rsp
         
             xmuf = scale
             xmur = xmuf
            xmuf2 = xmuf*xmuf 
c               AL = alphasPDF(xmur)
c               AL = 1.0d0
c              ALP = AL/2d0/Pi
              ALP = 1.0d0

            azmth = 2.0d0*pi
c               pf = 0.5d0*rsp
              ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth

c         write(*,*)'xa =', iplus, xa,x

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call setlum(f1,f2,xl)


c            call getPK(iplus,x,xmuf,p1,p2,p3,p4,SumPlus,SumReg,SumDel)
            call getPK(iplus,x,xmuf,p1,p2,p3,p4,SumPlus)
           

            if (iplus .eq. 1)  sig1 = xl(1)*SumPlus*coef
            if (iplus .eq. 0)  sig3 = xl(1)*SumPlus*coef
  
          if (iplus .eq. 1) then
            wgt1 = sig1/flux*ps2*xjac*vwgt
            PKplus_x = xnorm*wgt1/vwgt/2d0/eps

          elseif (iplus .eq.0) then
            wgt3 = sig3/flux*ps2*xjac*vwgt
            PKplus_1= xnorm*wgt3/vwgt/2d0/eps
          endif

           endif                ! bin choice
           enddo

            PK(k) = PKplus_x - PKplus_1

         enddo

        flo2_Plus = Alp * ( PK(1) + PK(2) )
      return
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
