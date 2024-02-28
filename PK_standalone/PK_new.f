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
      external Born_uU2eE
      call p2d_to_p1d_4(p,p1,p2,p3,p4)
        s12 = 2d0*dot(p1,p2)/x
        Cf = 4d0/3d0                      
        Tr = 0.5d0

      if (k .eq. 1) then   ! Cloice for [Leg1]       

        coefx = Al*(Cf*(-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {q,q,g}
        coefx = coefx*Born_uU2eE(0,p1,p2,p3,p4)

        coef1 = Al*(Cf*(-1)*Log(xmuf**2/(s12)))/2d0/Pi
        p1(0) = p1(0)/x
        p1(3) = p1(3)/x
        coef1 = coef1*Born_uU2eE(0,p1,p2,p3,p4) 

        ALLP(1) = Plusd((1+x**2)/(1-x))*(coefx - coef1)

      B=Al*Cf*((1d0+x)*Log((1d0-x)/x)+(1d0-x)-(1d0+x)*Log(1d0-x))/2d0/Pi
       Aplus=Plusd(2d0/(1d0-x)*Log((1d0-x)/x))+
     .    Plusd(2d0/(1d0-x)*dLog(1d0-x))
        call p2d_to_p1d_4(p,p1,p2,p3,p4)
        coefx = Al*Cf*Born_uU2eE(0,p1,p2,p3,p4)/2d0/Pi
        B = B * Born_uU2eE(0,p1,p2,p3,p4)
        p1(0) = p1(0)/x
        p1(3) = p1(3)/x
        coef1 = Al*Cf*Born_uU2eE(0,p1,p2,p3,p4)/2d0/Pi
        AllK(1)= Aplus*(coefx - coef1) - B 

      elseif (k .eq. 2 ) then  ! Cloice for [Leg2]

        coefx = Al*(Cf*(-1)*Log(xmuf**2/(s12*x)))/2d0/Pi    !{a,ai,i} => {q,q,g}
        coefx = coefx*Born_uU2eE(0,p1,p2,p3,p4)

        coef1 = Al*(Cf*(-1)*Log(xmuf**2/(s12)))/2d0/Pi
        p2(0) = p2(0)/x
        p2(3) = p2(3)/x
        coef1 = coef1*Born_uU2eE(0,p1,p2,p3,p4) 

        ALLP(1) = Plusd((1+x**2)/(1-x))*(coefx - coef1)



      B=Al*Cf*((1d0+x)*Log((1d0-x)/x)+(1d0-x)-(1d0+x)*Log(1d0-x))/2d0/Pi
       Aplus=Plusd(2d0/(1d0-x)*Log((1d0-x)/x))+
     .    Plusd(2d0/(1d0-x)*dLog(1d0-x))
        call p2d_to_p1d_4(p,p1,p2,p3,p4)
        coefx = Al*Cf*Born_uU2eE(0,p1,p2,p3,p4)/2d0/Pi
        B = B * Born_uU2eE(0,p1,p2,p3,p4)
        p2(0) = p2(0)/x
        p2(3) = p2(3)/x
        coef1 = Al*Cf*Born_uU2eE(0,p1,p2,p3,p4)/2d0/Pi
        AllK(1)= Aplus*(coefx - coef1) - B 

      endif
        SumP = 0d0 
        SumK = 0d0
       do i=1,1
        SumP = SumP + AllP(i)
        SumK = SumK + AllK(i)
       enddo
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
