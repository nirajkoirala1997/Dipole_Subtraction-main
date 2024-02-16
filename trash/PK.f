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
      call p2dtop1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)
        CLF = -1d0 ! Colour factor
        Al= 0.118d0
      if (k .eq. 1) then      ! Cloice for [Leg1]       
!       {a,ai,i} => {q,q,g}
        AllP(1) = Al*(Pqq(x)*CLF*dLog(xmuf**2/(s12*x)))/2d0/Pi
        AllP(2) = Al*(Pgq(x)*Clf*dLog(xmuf**2/(s12*x)))/2d0/Pi
        AllP(3) = Al*(Pgq(x)*Clf*dLog(xmuf**2/(s12*x)))/2d0/Pi !Contribution for anti-quark

        AllK(1) = Al*(Kbqq(x)-CLF*Kqq_til(x))/2d0/Pi !! #colour factor unclear
        AllK(2) = Al*(Kbgq(x)-CLF*Kgq_til(x))/2d0/Pi !! #colour factor unclear
        AllK(3) = Al*(Kbgq(x)-CLF*Kgq_til(x))/2d0/Pi !! #colour factor unclear

        SumP = 0d0 
        SumK = 0d0
       do i=1,3
        SumP = SumP + AllP(i)
        SumK = SumK + AllK(i)
       enddo

      elseif (k .eq. 2) then  ! Cloice for [Leg2]

        AllP(1) = Al*(Pqq(x)*CLF*dLog(xmuf**2/(s12*x)))/2d0/Pi
        AllP(2) = Al*(Pgq(x)*Clf*dLog(xmuf**2/(s12*x)))/2d0/Pi
        AllP(3) = Al*(Pgq(x)*Clf*dLog(xmuf**2/(s12*x)))/2d0/Pi !Contribution for anti-quark

        AllK(1) = Al*(Kbqq(x)-CLF*Kqq_til(x))/2d0/Pi !! #colour factor unclear
        AllK(2) = Al*(Kbgq(x)-CLF*Kgq_til(x))/2d0/Pi !! #colour factor unclear
        AllK(3) = Al*(Kbgq(x)-CLF*Kgq_til(x))/2d0/Pi !! #colour factor unclear

        SumP = 0d0 
        SumK = 0d0
       do i=1,3
        SumP = SumP + AllP(i)
        SumK = SumK + AllK(i)
       enddo
      endif

      return
      end
c-----------------------------------------------------------  !P-terms
      function Pqg(x) ! OR Pqbg antiquark also same                                    
      implicit double precision (a-h,o-z)
        Cf = 4d0/3d0
         Pqg = Cf*(1d0+(1d0-x)**2)/x
      return
      end

      function Pgq(x)
      implicit double precision (a-h,o-z)
         Tr = 0.5d0
         Pgq = Tr*(x**2+(1-x)**2)
      return
      end

      function Pqq(x)
      implicit double precision (a-h,o-z)
        Cf = 4d0/3d0
        Pqq = Cf*Plusd((1+x**2)/(1-x))
      return
      end

      function Pgg(x)
      implicit double precision (a-h,o-z)
         Ca = 3d0
         Nf = 1d0
         Tr = 0.5d0
         Pgg = 2d0*Ca*(Plusd(1d0/(1d0-x)) + (1d0-x)/x -1d0+ x*(1d0-x)) +
     .         delta(1d0-x)*(11d0*Ca/6d0 - 2d0*Nf*Tr/3d0)
      return
      end
c-----------------------------------------------------------!K-terms

      function Kbqq(x)
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
        Cf = 4d0/3d0
        Kbqq = Cf*(Plusd(2d0/(1d0-x)*dLog((1d0-x)/x))-
     .        (1d0-x)*dLog((1d0-x)/x)+(1d0-x)-
     .         delta(1d0-x)*(5d0-Pi**2))
      return
      end                

      function Kbgq(x) !Kb[GQ]
      implicit double precision (a-h,o-z)
        Tr = 0.5d0 
        Kbgq = Pgq(x)*dLog((1d0-x)/x) + Tr*2d0*x*(1d0-x)
      return
      end                

      function Kbqg(x) !Kb[QG]
      implicit double precision (a-h,o-z)
        Tr = 0.5d0 
        Cf = 4d0/3d0
        Kbqg = Pqg(x)*dLog((1d0-x)/x) + Cf*x 
      return
      end                

      function Kqq_til(x)
      implicit double precision (a-h,o-z)
        Cf = 4d0/3d0
        Kqq_til = Pqq(x)*dLog(1d0-x) + 
     .    Cf*(Plusd(2d0/(1d0-x)*dLog(1d0-x))-Pi**3*delta(1d0-x)/3d0)
      return
      end

      function Kgq_til(x)
      implicit double precision (a-h,o-z)
        Kgq_til = Pgq(x)*dLog(1d0-x)  
      return
      end
c----------------------------------------------------------- !Misc-functions
                                                            
      function Plusd(xx)
      implicit double precision (a-h,o-z)
        Plusd = xx
      return
      end

      function delta(xx)
      implicit double precision (a-h,o-z)
        delta = 0d0
      return
      end
