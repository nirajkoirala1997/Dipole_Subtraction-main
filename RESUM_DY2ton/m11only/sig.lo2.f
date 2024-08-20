c----------------------------------------C
c     2 -> 2 lo partonic cross section  -C
c----------------------------------------C

      subroutine sig_lo2(p1,p2,p3,p4,sig)
      implicit double precision (a-h,o-z)
      dimension f1(-6:6),f2(-6:6)
      dimension xl(15)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension sig(10)
c      dimension sigma(10)
      common/stu/s,t,u
      COMMON/FLAGS/IFLAG
      COMMON/VMASS/XMZ,XMW,XMH,GMZ,GMW

      s=2d0*dot(p1,p2)
      t=xmz**2-2d0*dot(p1,p3)
      u=-(s+t)+2.0d0*xmz**2

      call setsigma(s,t,u,sig)
      return
      end

