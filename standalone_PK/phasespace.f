cc +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
cc Two body phase space
cc +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
c
c      subroutine kinvar2(xx,xxinvmass,p1,p2,p3,p4)
c      implicit double precision (a-h,o-z)
c      dimension xx(4)
c      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
c      dimension qa(0:3),qb(0:3)
c      common/angle/cststar
c
c      s=1000d0**2
c      v=xx(1)
c      omv=1d0-v
c      srs2=0.5*sqrt(s)
c
cc     incoming parton 4-vectors
c      p1(0)=srs2
c      p1(1)=0d0
c      p1(2)=0d0
c      p1(3)=p1(0)
c
c      p2(0)=srs2
c      p2(1)=0d0
c      p2(2)=0d0
c      p2(3)=-p2(0)
c
cc     outgoing parton 4-vectors
c      p3(0)=srs2*(v+omv)
c      p3(1)=sqrt(s*v*omv)
c      p3(2)=0d0
c      p3(3)=srs2*(v-omv)
c
c      p4(0)=p1(0)+p2(0)-p3(0)
c      p4(1)=p1(1)+p2(1)-p3(1)
c      p4(2)=p1(2)+p2(2)-p3(2)
c      p4(3)=p1(3)+p2(3)-p3(3)
c
cc     p3 + p4
c      q(0) = p3(0) +p4(0)
c      q(1) = p3(1) +p4(1)
c      q(2) = p3(2) +p4(2)
c      q(3) = p3(3) +p4(3)
c
cc     invariant mass of diphoton pair.
c      s34    = 2*dot(p3,p4)
c      xxinvmass= sqrt(s34)
c
cc     cos-theta
c      ccst1=p3(3)/p3(0)
c      ccst2=p4(3)/p4(0)
c
cc     ppt1 ppt2
c      ppt1=sqrt(p3(1)**2+p3(2)**2)
c      ppt2=sqrt(p4(1)**2+p4(2)**2)
c
c      qa(0) =  p3(0) - p4(0)
c      qa(1) =  p3(1) - p4(1)
c      qa(2) =  p3(2) - p4(2)
c      qa(3) =  p3(3) - p4(3)
c
c      qb(0) =  p3(0) + p4(0)
c      qb(1) =  p3(1) + p4(1)
c      qb(2) =  p3(2) + p4(2)
c      qb(3) =  p3(3) + p4(3)
c
c      p1qa = dot(p1,qa)
c      p1qb = dot(p1,qb)
c
c      cststar = p1qa/p1qb
c      return
c      end
      subroutine kinvar2(xx,xxinvmass,p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension xx(10)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3)
      dimension qa(0:3),qb(0:3)
      common/energy/s

      xa=xx(1)
      xb=xx(2)
      v=xx(3)
      omv=1d0-v     
      
c      s=s*xa*xb 
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
      p3(0)=srs2*(xa*v+xb*omv)
      p3(1)=dsqrt(s*xa*xb*v*omv)
      p3(2)=0d0
      p3(3)=srs2*(xa*v-xb*omv)

      p4(0)=p1(0)+p2(0)-p3(0)
      p4(1)=p1(1)+p2(1)-p3(1)
      p4(2)=p1(2)+p2(2)-p3(2)
      p4(3)=p1(3)+p2(3)-p3(3)

c     p3 + p4
      q(0) = p3(0) +p4(0)
      q(1) = p3(1) +p4(1)
      q(2) = p3(2) +p4(2)
      q(3) = p3(3) +p4(3)

c     invariant mass of diphoton pair.
      s34    = 2*dot(p3,p4)
      xxinvmass= dsqrt(s34)

c     yy1: rapidity of photons 3 
      yrpda  = (p3(0)+p3(3))/(p3(0)-p3(3))
      yy1    = 0.5*dlog(yrpda)

c     yy2: rapiditiy of photon 4
      yrpdb  = (p4(0)+p4(3))/(p4(0)-p4(3))
      yy2    = 0.5*dlog(yrpdb)
      
c     rapidity YY
      p2q = dot(p2,q)
      p1q = dot(p1,q)
      rr  = xa*p2q/(xb*p1q)
      YY12  = dlog(rr)/2.d0

c     cos-theta
      ccst1=p3(3)/p3(0)
      ccst2=p4(3)/p4(0)

c     ppt1 ppt2
      ppt1=dsqrt(p3(1)**2+p3(2)**2)
      ppt2=dsqrt(p4(1)**2+p4(2)**2)

      qa(0) =  p3(0) - p4(0)
      qa(1) =  p3(1) - p4(1)
      qa(2) =  p3(2) - p4(2)
      qa(3) =  p3(3) - p4(3)

      qb(0) =  p3(0) + p4(0)
      qb(1) =  p3(1) + p4(1)
      qb(2) =  p3(2) + p4(2)
      qb(3) =  p3(3) + p4(3)

      p1qa = dot(p1,qa)
      p1qb = dot(p1,qb)

      cststar = p1qa/p1qb

      return
      end
c---------------------------------------------------------------------
         subroutine p2d_to_p1d_4(p,p1,p2,p3,p4)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
        do i=0,3
         p1(i)=p(i,1)
         p2(i)=p(i,2)
         p3(i)=p(i,3)
         p4(i)=p(i,4)
        enddo
         end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
         subroutine p1d_to_p2d_4(p1,p2,p3,p4,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
         p(i,4)=p4(i)
        enddo
         end
c---------------------------------------------------------------------
 
