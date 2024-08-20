      subroutine event(x1,x2,x4,flo2,virt,virtN)
c      function flo2(x1,x2,x4)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension sig(10)
      dimension f1(-6:6),f2(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.38937966d12)          ! cross section will be in fb
      
      COMMON/MACHINE/RS
c      common/momenta/p1,p2,p3,p4
c      common/bin/xq,xeps,xlow,xhigh
      common/factscale/xmuf
      common/param/aem,xmur,lambda
      common/xmcoeff/xc1,xc2
      common/isub/io,is
      common/angle/cststar
      common/amf/am3,am4,am5
c      common/mass/z,h
      COMMON/Qvalue/Qmin,Qmax,XEPS
      COMMON/VMASS/XMZ,XMW,XMH,GMZ,GMW

c      common/scale/scale,zz
      
      xa     = x1
      xb     = x2

c      z = 91.1876d0**2
c      h = 0d0
      am3 = xmz
      am4 = xmz
c      am3 = dsqrt(z) 
c      am4 = dsqrt(z)

      call kinvar2(x1,x2,x4,xinvmass,pf,p1,p2,p3,p4)
c      call kinvar2(x1,x2,x4,xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,pf)
c      call cuts0(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,ipass)

      s = RS*RS
      xinvmass2 = 2.0d0*xmz**2 + 2.0d0*dot(p3,p4)
      scale  = dsqrt(xinvmass2)

c      TAUH = Scale*scale/s
c      zz = TAUH/X1/X2
c      scale2  = scale*scale

      if (scale.ge.(am3+am4)) then
      ipass = 1
      else
      ipass = 0
      endif

      IF(ipass.eq.1)THEN

            call sig_lo2(p1,p2,p3,p4,sig)

            pfpi = 2d0*pf/dsqrt(s*xa*xb)
            xnorm=hbarc2*pfpi/16d0/pi/(s*xa*xb)
            wgt=xnorm*sig(1)
            flo2=wgt
            virt = sig(2)
            virtN = sig(4)

            return              ! There is a factor of 1/2
      ELSEIF(ipass.eq.0)THEN
         flo2=0d0
         virt=0d0
         virtN=0d0
         return
      ENDIF

      end
