c Generated by AutoDipole
c Kouhei Hasegawa, Sven Moch, and Peter Uwer, 2009
c Filename:Iterm.f
c Subroutine: Iterm evaluates all I-terms
c
c input
c   real p      :a phase space point
c output
c   real coef   :coefficient of all I-terms
c   real SumI   :sum of all coefficients
c 
c Subroutine to calculate I terms 
c LO and Virtual process:{{u, p[1]}, {ubar, p[2]}} --> {{e, p[3]}, {ebar, p[4]}}
c I(i) = coef(i,-2)*eps^(-2) + coef(i,-1)*eps^(-1) + coef(i,0)  
c SumI[-2,-1,0] = Sum_{i=1}^{2}coef[-2,-1,0] 

      subroutine Iterm(p,coef,SumI) 
      implicit none 

      integer i,j 
      double precision p(0:3,1:4),coef(2,-2:0),SumI(-2:0)
      double precision Pi,rtwo,Eul,AL,CF,CA,TR,mt,mb,dot,mu,Nf,log2
      double precision CLV(2),q(0:3,1:4),dilog
      double precision s12,s13,s14,s23,s24,s34 
 
      common /MASS/ mt,mb  
      common /usedalpha/ AL  
      common/renor_scale/mu
    
      Pi=3.141592653589793238D0 
      Eul=0.5772156649015328606065120d0 
      rtwo=dsqrt(2.d0) 
      log2=Log(2.d0) 
      CF=4.D0/3.D0 
      CA=3.D0 
      TR=0.5D0 
c      mu=174.3D0 
      Nf=1D0 
      
      s12=2.d0*dot(p(0,1),p(0,2)) 
      s13=2.d0*dot(p(0,1),p(0,3)) 
      s14=2.d0*dot(p(0,1),p(0,4)) 
      s23=2.d0*dot(p(0,2),p(0,3)) 
      s24=2.d0*dot(p(0,2),p(0,4)) 
      s34=2.d0*dot(p(0,3),p(0,4)) 
  
      do i=0,3 
       q(i,1)=p(i,1) 
       q(i,2)=p(i,2) 
      enddo 
      
      call cmatrix1(p(0,1),1,CLV(1)) 
      
      call cmatrix1(p(0,1),2,CLV(2)) 
      
      coef(1,-2)=
     - -0.5*(AL*CLV(1))/Pi
   
      coef(1,-1)=
     -         (AL*CLV(1)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s12)))/(4.*Pi)
   
      coef(1,0)=
     -         (AL*CLV(1)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s12) - 
     -      6*Log(mu**2/s12)**2))/(24.*Pi)
   
      coef(2,-2)=
     - -0.5*(AL*CLV(2))/Pi
   
      coef(2,-1)=
     -         (AL*CLV(2)*(-3 + 2*Eul - 2*Log(4*Pi) - 
     -      2*Log(mu**2/s12)))/(4.*Pi)
   
      coef(2,0)=
     -         (AL*CLV(2)*(-60 - 6*Eul**2 + 7*Pi**2 - 18*Log(4*Pi) - 
     -      6*Log(4*Pi)**2 + 6*Eul*(3 + 2*Log(4*Pi)) + 
     -      6*(-3 + 2*Eul - 2*Log(4*Pi))*Log(mu**2/s12) - 
     -      6*Log(mu**2/s12)**2))/(24.*Pi)
   
   
      SumI(-2) = 0.d0 
      SumI(-1) = 0.d0 
      SumI(-0) = 0.d0 
 
      do i=1,2 
      SumI(-2) = SumI(-2) + coef(i,-2) 
      SumI(-1) = SumI(-1) + coef(i,-1) 
      SumI(0) = SumI(0) + coef(i,0) 
      enddo 
 
 
      return 
      end 

      subroutine mexchange(q,p,i,j)
      implicit none
      integer i,j,k,l
      double precision p(0:3,1:4),q(0:3,1:4)
      do k=0,3
       do l=3,4
        q(k,l) = p(k,l)
       enddo
       q(k,i) = p(k,j)
       q(k,j) = p(k,i)
      enddo
      return
      end
 
      subroutine mmap(q,p)
      implicit none
      integer i
      double precision  p(0:3),q(0:3)
      do i=0,3
       q(i) = p(i)
      enddo
      return
      end
 
 
