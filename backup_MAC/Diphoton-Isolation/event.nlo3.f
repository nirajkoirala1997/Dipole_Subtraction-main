      function fnlo3(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
      dimension f1(-6:6),f2(-6:6)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=0.3894d9)

      common/energy/s
      common/momenta5/p1,p2,p3,p4,p5
      common/bin/xq,xeps,xlow,xhigh
      common/cfbin/cf,cfeps,clow,chigh
      common/yfbin/yf,yfeps,ylow,yhigh
      common/xrange/xrlow,xrhigh
      common/factscale/xmuf
      common/param/aem,xmur,lambda
      common/xmcoeff/xc1,xc2
      common/isub/io,is
      common/angle/cststar
 

c        yy(1)=   4.1782791168462528E-002
c        yy(2)=   1.7352003311930358E-002
c        yy(3)=  0.28579477634085454     
c        yy(4)=  0.47774661455815909     
c        yy(5)=  0.67105779080538752     
c        yy(6)=  0.99465843471179438     





      call kinvar3(yy,xjac,xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,
     &     r34,r35,r45,Et5,QT)
      
      call cuts3(xinvmass,y1,y2,Y,cst1,cst2,pt1,pt2,
     &     r34,r35,r45,Et5,QT,ipass)
c      
c	print*,"P1:",p1
c	print*,"P2:",p2
c	print*,"P3:",p3
c	print*,"P4:",p4
c	print*,"P5:",p5
c	stop




      xa =yy(1)
      xb =yy(2)
      v  =yy(3)

	rsp = dsqrt(xa*xb*s)
      scale = xinvmass
      pobl = cststar

      IF(ipass.eq.1)THEN

         if(   scale.gt.xlow .and. scale.lt.xhigh
c         if(   scale.gt.xrlow .and. scale.lt.xrhigh
c     &  .and. pobl.ge.clow .and. pobl.le.chigh
     &        )then

c        write(*,*)'Qmin, Qmax', xlow, xhigh
c            write(*,*)'io,is',io,is
c            write(*,*)'Q =', xinvmass

            xmuf=xc1*scale
            xmur=xc2*scale

            call pdf(xa,xmuf,f1)
            call pdf(xb,xmuf,f2)
            call sig_nlo3(f1,f2,p1,p2,p3,p4,p5,sig)
c	print*,yy
c	print*,"sigma:",sig
c	stop

            xnorm=hbarc2*v/64/4/4/(pi**4)
            wgt=xjac*xnorm*sig*vwgt
            fnlo3=wgt/vwgt/2.d0/xeps

c          pi_1 = 0.5d0*rsp
c          flux = 4d0*pi_1*rsp
c          xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
c          wgt=xjac*xnorm*sig*vwgt
c          fnlo3=wgt/vwgt/2d0/xeps


c       print*,"xx(1)=",yy(1)
c       print*,"xx(2)=",yy(2)
c       print*,"xx(3)=",yy(3)
c       print*,"xx(4)=",yy(4)
c       print*,"xx(5)=",yy(5)
c       print*,"xx(6)=",yy(6)

c   	print*,"Integrand:",fnlo3
c	print*," "
c	stop
            return              !There is a factor of 1/2
         else                   !from xeps interval
            fnlo3=0d0
            return
         endif         
      ELSEIF(ipass.eq.0)THEN
         fnlo3=0d0
         return
      ENDIF

      end


