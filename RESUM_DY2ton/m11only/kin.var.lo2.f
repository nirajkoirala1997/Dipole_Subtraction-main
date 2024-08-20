c      subroutine kinvar2(x1,x2,x4,
c     &     xxinvmass,yy1,yy2,YY12,ccst1,ccst2,ppt1,ppt2,pf) 
      subroutine kinvar2(x1,x2,x4,xinvmass,pf,p1,p2,p3,p4)
c      subroutine kinvar2(x1,x2,x4,p1,p2,p3,p4) 
      implicit double precision (a-h,o-z)
      dimension xx(10)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      COMMON/MACHINE/RS
c      common/momenta/p1,p2,p3,p4
      common/add_par/xms,nd
      common/amf/am3,am4,am5
      common/angle/cststar

      xa=x1
      xb=x2
      v=x4
      omv=1d0-v     

      s = RS*RS
      sp = xa*xb*s
c      write(*,*)'s,sp = ',xa,xb,sp
      rsp = dsqrt(sp)

      srs2=0.5*dsqrt(s)

c     incoming parton 4-vectors
      p1(0)=srs2*xa
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=p1(0)

      p2(0)=srs2*xb
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=-p2(0)

c     outgoing parton 4-vectors
c     For massless case
c      p3(0)=srs2*(xa*v+xb*omv)
c      p3(1)=dsqrt(s*xa*xb*v*omv)
c      p3(2)=0d0
c      p3(3)=srs2*(xa*v-xb*omv)
c      write(*,*)'m3, m4 =', am3, am4

c For massive case where m3 and m4 are non-zero
c ---------------------------------------------
c       Parametrization in the c.o.m. frame
c ---------------------------------------------
        xmd = (am3**2 - am4**2)
        e3 = 0.5d0*(rsp + xmd/rsp)
        pf = dsqrt(e3**2 - am3**2)
        ct = 2d0*v-1
        st = dsqrt(1d0 - ct**2d0)
        px3 = pf*st
        py3 = 0d0
        pz3 = pf*ct
        
c ---------------------------------------------
c       Boost to lab frame
c ---------------------------------------------
        beta=(xb-xa)/(xa+xb)
        gamma=1d0/dsqrt(1d0-beta*beta)

        p3(0) = gamma*(e3 - beta*pz3)
        p3(1) = px3
        p3(2) = py3
        p3(3) = gamma*(pz3 - beta*e3)

      p4(0)=p1(0)+p2(0)-p3(0)
      p4(1)=p1(1)+p2(1)-p3(1)
      p4(2)=p1(2)+p2(2)-p3(2)
      p4(3)=p1(3)+p2(3)-p3(3)

      return
      end


