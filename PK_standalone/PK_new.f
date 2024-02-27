! a -> ai+ i ||  q -> q g || P{a,ai} -> P{q,q} || a- Incoming
!                                               | ai- Undergoing Born
!                                               |    process               
! Kb = Kbar      
c----------------------------------------------------------- 
      subroutine getPK(k,x,xmuf,p,SumP,SumK)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      dimension AllP(5),AllK(5)
      dimension  p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
      common /usedalpha/ AL      
      call p2d_to_p1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
        Cf = 4d0/3d0                      
        Tr = 0.5d0

      if (k .eq. 1 .or. k .eq. 3) then   ! Cloice for [Leg1]       
!       
        AllP(1) = Al*(Cf*Plusd((1+x**2)/(1-x))*
     .            (-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {q,q,g}
        AllP(2) = Al*(Tr*(2*x**2+1-
     .                 2*x)*(-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {g,q,q~}

        AllK(1) = Al*(Cf*(Plusd(2d0/(1d0-x)*Log((1d0-x)/x))-
     .        (1d0+x)*Log((1d0-x)/x)+(1d0-x)-
     .         delta(1d0-x)*(5d0-Pi**2))- Cf*(1d0+x)*Log(1d0-x) + 
     .    Cf*(Plusd(2d0/(1d0-x)*dLog(1d0-x))-
     .              Pi**2*delta(1d0-x)/3d0))/2d0/Pi 

        AllK(2) = Al*(Tr*(x**2+(1-x)**2)*Log((1d0-x)/x) +
     .           Tr*2d0*x*(1d0-x)-(-1)*Tr*(x**2+
     .              (1-x)**2)*Log(1d0-x))/2d0/Pi  
        

      elseif (k .eq. 2 .or. k .eq. 3) then  ! Cloice for [Leg2]


        AllP(1) = Al*(Cf*Plusd((1+x**2)/(1-x))*
     .            (-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {q,q,g}
        AllP(2) = Al*(Tr*(2*x**2+1-
     .                 2*x)*(-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {g,q,q~}

        AllK(1) = Al*(Cf*(Plusd(2d0/(1d0-x)*Log((1d0-x)/x))-
     .        (1d0+x)*Log((1d0-x)/x)+(1d0-x)-
     .         delta(1d0-x)*(5d0-Pi**2))- Cf*(1d0+x)*Log(1d0-x) + 
     .    Cf*(Plusd(2d0/(1d0-x)*dLog(1d0-x))-
     .              Pi**2*delta(1d0-x)/3d0))/2d0/Pi !! #colour factor unclear

        AllK(2) = Al*(Tr*(x**2+(1-x)**2)*Log((1d0-x)/x) +
     .           Tr*2d0*x*(1d0-x)-(-1)*Tr*(x**2+
     .              (1-x)**2)*Log(1d0-x))/2d0/Pi  
        
      endif

c        if (k .eq. 3) then
c                print*,"gaya"
c                do i=1,2
c                print*,Allp(i),allk(i)
c                print*,"LOG:",Log(xmuf**2/(s12*x))
c                print*,"next:",xmuf**2,s12*x,x
c                enddo
c        endif
        SumP = 0d0 
        SumK = 0d0
       do i=1,2
        SumP = SumP + AllP(i)
        SumK = SumK + AllK(i)
       enddo
       if (x .eq. 1d0) then
               sump=0d0
               sumk=0d0
       endif
c       print*,"P:",SumP
c       print*,"K:",SumK

      return
      end
c-----------------------------------------------------------  !P-terms
c----------------------------------------------------------- !Misc-functions
c                                                            
c      function Plusd(xx)
c      implicit double precision (a-h,o-z)
c        Plusd = xx
c      return
c      end
c
c      function delta(xx)
c      implicit double precision (a-h,o-z)
c        delta = 0d0
c      return
c      end
