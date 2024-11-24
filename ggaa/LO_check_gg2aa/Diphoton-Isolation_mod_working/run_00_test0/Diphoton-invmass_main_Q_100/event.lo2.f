      function flo2(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)
      
      common/energy/s
      common/momenta/p1,p2,p3,p4
      common/bin/xq,xeps,xlow,xhigh
      common/cfbin/cf,cfeps,clow,chigh
      common/yfbin/yf,yfeps,ylow,yhigh
      common/xrange/xrlow,xrhigh
      common/factscale/xmuf
      common/param/aem,xmur,lambda
      common/xmcoeff/xc1,xc2
      common/isub/io,is
      common/angle/cststar

      tau = xq*xq/s

      xamin = tau
      xamax = 1.0d0
      xajac = (xamax - xamin)
      xa     = xamin+yy(1)*xajac
      xb = tau/xa
      xc = yy(2)

c      xa = yy(1)
c      xb = yy(2)
c      xc = yy(3)

      call kinvar2(xa,xb,xc,xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2)
      call cuts0(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,ipass)
c      ipass = 1

      scale  = xinvmass

      pobl = cststar

      IF(ipass.eq.1)THEN

c         if(   scale.gt.xlow .and. scale.lt.xhigh
c         if(   scale.gt.xrlow .and. scale.lt.xrhigh
c     &  .and. pobl.ge.clow .and. pobl.le.chigh
c     &        )then

c           write(*,*)scale, pobl, cf
c            xmuf=xc1*scale
            xmuf=xc1*scale
            xmur=xc2*scale

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call sig_lo2(f1,f2,p1,p2,p3,p4,sig)

            xnorm=hbarc2/16d0/pi/(s*xa*xb)
            wgt=xnorm*sig*vwgt
c            flo2=wgt/vwgt/2.d0/xeps
            flo2=wgt/vwgt*xajac*2*xq/s/xa
            return              ! There is a factor of 1/2
c         else                   ! from xeps interval
c            flo2=0d0
c            return
c         endif
      ELSEIF(ipass.eq.0)THEN
         flo2=0d0
         return
      ENDIF

      end
