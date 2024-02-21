c Generated by AutoDipole
c Kouhei Hasegawa, Sven Moch, and Peter Uwer, 2009
c Filename:PKterm.f
c Subroutine: PKterm evaluates all P and K terms
c
c input
c   real p      :a phase space point
c   real x      :a momentum fraction
c output
c   real SumP   :sum of all P-terms
c   real SumK   :sum of all K-terms

c Subroutine to calculate P and K terms of Leg 1
c LO and Virtual process:{{u, p[1]}, {ubar, p[2]}} --> {{e, p[3]}, {ebar, p[4]}}
c B1:1-2
c B4u:3-4
c SumP = Sum_{i=1}^{4}P(i) 
c SumK = Sum_{i=1}^{4}K(i) 

      subroutine PKterm1(p,x,SumP,SumK) 
      implicit none 

      integer i,j 
      double precision p(0:3,1:4),SumP,SumK,Pt(4),K(4)
      double precision SumP1,SumP3u,SumP3ubar,SumP4u,SumP4ubar
      double precision SumK1,SumK3u,SumK3ubar,SumK4u,SumK4ubar
      double precision Pi,rtwo,Eul,AL,CF,CA,TR,mt,mb,dot,mu,Nf,log2
      double precision muF,x,plusd,delta
      double precision CLPK(4),q(0:3,1:4),dilog
      double precision s12,s13,s14,s23,s24,s34 
 
      common /MASS/ mt,mb  
      common /usedalpha/ AL  
      common/factscale/muf
    
      Pi=3.141592653589793238D0 
      Eul=0.5772156649015328606065120d0 
      rtwo=dsqrt(2.d0) 
      log2=Log(2.d0) 
      CF=4.D0/3.D0 
      CA=3.D0 
      TR=0.5D0 
      Nf=1.D0 


      s12=2.d0*dot(p(0,1),p(0,2))/x 
      s13=2.d0*dot(p(0,1),p(0,3))/x 
      s14=2.d0*dot(p(0,1),p(0,4))/x 
      s23=2.d0*dot(p(0,2),p(0,3)) 
      s24=2.d0*dot(p(0,2),p(0,4)) 
      s34=2.d0*dot(p(0,3),p(0,4)) 
  
      do i=0,3 
       q(i,1)=p(i,1) 
       q(i,2)=p(i,2) 
      enddo 
      
      call smatrix(p(0,1),CLPK(1)) 
      call cmatrix1(p(0,1),1,CLPK(2)) 

      call smatrix(p(0,1),CLPK(3)) 
      call cmatrix1(p(0,1),1,CLPK(4)) 

      Pt(1)=
     - 0
   
      Pt(2)=
     -         (AL*CLPK(2)*Log(muF**2/(s12*x))*
     -    plusd((1 + x**2)/(1 - x)))/(2.*Pi)
   
      Pt(3)=
     - 0
   
      Pt(4)=
     -         (AL*TR*((1 - x)**2 + x**2)*CLPK(4)*Log(muF**2/(s12*x)))/
     -  (2.*CF*Pi)
   


      K(1)=
     -         (AL*CLPK(1)*(-(CF*(5 - Pi**2)*delta(1 - x)) + 
     -      CF*(1 - x - (1 + x)*Log((1 - x)/x) + 
     -         plusd((2*Log((1 - x)/x))/(1 - x)))))/(2.*Pi)
   
      K(2)=
     -         -0.5*(AL*CLPK(2)*(-(CF*(1 + x)*log(1 - x)) + 
     -       CF*(-0.3333333333333333*(Pi**2*delta(1 - x)) + 
     -          plusd((2*Log(1 - x))/(1 - x)))))/(CF*Pi)
   
      K(3)=
     -         (AL*CLPK(3)*(2*TR*(1 - x)*x + 
     -      TR*((1 - x)**2 + x**2)*Log((1 - x)/x)))/(2.*Pi)
   
      K(4)=
     - -0.5*(AL*TR*(1 - 2*x + 2*x**2)*CLPK(4)*Log(1 - x))/(CF*Pi)
   
      SumP1 = 0.d0 
      SumK1 = 0.d0 
      do i=1,2 
      SumP1 = SumP1 + Pt(i) 
      SumK1 = SumK1 + K(i) 
      enddo

      SumP4u = 0.d0 
      SumK4u = 0.d0 
      do i=3,4 
      SumP4u = SumP4u + Pt(i) 
      SumK4u = SumK4u + K(i) 
      enddo

      SumP = 0.d0 
      SumP = SumP1+SumP4u
      SumK = 0.d0 
      SumK = SumK1+SumK4u
      if ((SumP .ne. SumP) .or. (SumK .ne. SumK)) then
        SumP = 0d0
        SumK = 0d0
      endif
      return 
      end 

c Subroutine to calculate P and K terms of Leg 2
c LO and Virtual process:{{u, p[1]}, {ubar, p[2]}} --> {{e, p[3]}, {ebar, p[4]}}
c B1:1-2
c B4ubar:3-4
c SumP = Sum_{i=1}^{4}P(i) 
c SumK = Sum_{i=1}^{4}K(i) 

      subroutine PKterm2(p,x,SumP,SumK) 
      implicit none 

      integer i,j 
      double precision p(0:3,1:4),SumP,SumK,Pt(4),K(4)
      double precision SumP1,SumP3u,SumP3ubar,SumP4u,SumP4ubar
      double precision SumK1,SumK3u,SumK3ubar,SumK4u,SumK4ubar
      double precision Pi,rtwo,Eul,AL,CF,CA,TR,mt,mb,dot,mu,Nf,log2
      double precision muF,x,plusd,delta
      double precision CLPK(4),q(0:3,1:4),dilog
      double precision s12,s13,s14,s23,s24,s34 
 
      common /MASS/ mt,mb  
      common /usedalpha/ AL  
      common/factscale/muf
    
      Pi=3.141592653589793238D0 
      Eul=0.5772156649015328606065120d0 
      rtwo=dsqrt(2.d0) 
      log2=Log(2.d0) 
      CF=4.D0/3.D0 
      CA=3.D0 
      TR=0.5D0 
      Nf=1.D0 
      AL=0.118d0


      s12=2.d0*dot(p(0,1),p(0,2))/x 
      s13=2.d0*dot(p(0,1),p(0,3)) 
      s14=2.d0*dot(p(0,1),p(0,4)) 
      s23=2.d0*dot(p(0,2),p(0,3))/x 
      s24=2.d0*dot(p(0,2),p(0,4))/x 
      s34=2.d0*dot(p(0,3),p(0,4)) 
  
      do i=0,3 
       q(i,1)=p(i,1) 
       q(i,2)=p(i,2) 
      enddo 
      
      call smatrix(p(0,1),CLPK(1)) 
      call cmatrix1(p(0,1),2,CLPK(2)) 

      call smatrix(p(0,1),CLPK(3)) 
      call cmatrix1(p(0,1),2,CLPK(4)) 
      Pt(1)=
     - 0
   
      Pt(2)=
     -         (AL*CLPK(2)*Log(muF**2/(s12*x))*
     -    plusd((1 + x**2)/(1 - x)))/(2.*Pi)
   
      Pt(3)=
     - 0
   
      Pt(4)=
     -         (AL*TR*((1 - x)**2 + x**2)*CLPK(4)*Log(muF**2/(s12*x)))/
     -  (2.*CF*Pi)
   


      K(1)=    ! this is Kbar_qq term
     -         (AL*CLPK(1)*(-(CF*(5 - Pi**2)*delta(1 - x)) + 
     -      CF*(1 - x - (1 + x)*Log((1 - x)/x) + 
     -         plusd((2*Log((1 - x)/x))/(1 - x)))))/(2.*Pi)
   
      K(2)=    !this is the Ktil_qq term
     -         -0.5*(AL*CLPK(2)*(-(CF*(1 + x)*log(1 - x)) + 
     -       CF*(-0.3333333333333333*(Pi**2*delta(1 - x)) + 
     -          plusd((2*Log(1 - x))/(1 - x)))))/(CF*Pi)
   
      K(3)=    ! this is K_gq term
     -         (AL*CLPK(3)*(2*TR*(1 - x)*x + 
     -      TR*((1 - x)**2 + x**2)*Log((1 - x)/x)))/(2.*Pi)
   
      K(4)=    ! this is P+regular term in Ktil_gq
     - -0.5*(AL*TR*(1 - 2*x + 2*x**2)*CLPK(4)*Log(1 - x))/(CF*Pi)
   
      SumP1 = 0.d0 
      SumK1 = 0.d0 
      do i=1,2 
      SumP1 = SumP1 + Pt(i) 
      SumK1 = SumK1 + K(i) 
      enddo

      SumP4ubar = 0.d0 
      SumK4ubar = 0.d0 
      do i=3,4 
      SumP4ubar = SumP4ubar + Pt(i) 
      SumK4ubar = SumK4ubar + K(i) 
      enddo

      SumP = 0.d0 
      SumP = SumP1+SumP4ubar
      SumK = 0.d0 
      SumK = SumK1+SumK4ubar
      if ((SumP .ne. SumP) .or. (SumK .ne. SumK)) then
        SumP = 0d0
        SumK = 0d0
      endif

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
 
      double precision function plusd(xx)
      implicit none
      double precision xx
      plusd=xx
      return
      end

      double precision function delta(yy)
      implicit none
      double precision yy
      delta=0.d0
      return
      end

