      function alphas(nkind,nord,xmz,pdfname,alphsmz,mur)

      implicit none
      integer nkind, nord,mem,i,nmax
      real*8 as,mur,xmz,alphsmz,as_mz,as_n3loxs,step,PI
      real*8 alphasPDF,asopi,alphas
      character pdfname*100
      external as_n3loxs
      DATA PI/3.141592653589793238462643D0/


      if (nkind .eq. 1)  then             ! LHAPDF alphas
      alphas = alphasPDF(mur)

      elseif (nkind .eq. 2)  then             ! n3loxs alphas
      as_mz = alphsmz/Pi
      asopi = as_n3loxs(mur,nord,xmz,as_mz)
      alphas = asopi*Pi
c      write(*,*)'alphas =', mur,nord,xmz,xmz,as_mz,alphas
      endif

      return
      end


       
        function as_n3loxs(mur, N, MZ, as_mz)
        implicit none
        integer nsteps, N, i
        real*8 MZ,Nc,nf,Zeta3,Zeta5,b0,b1,b2,b3
        real*8 steps,res,mur,as_mz
        real*8  k1, k2, k3, k4
        real*8 as_n3loxs

c       QCD parameters:
        Nc = 3.0d0
        nf = 5.0d0

c       Zeta function
        Zeta3 = 1.2020569031595942853997381615114499907649862923405d0
        Zeta5 = 1.0369277551433699263313654864570341680570809195019d0

c       QCD beta function
        b0 = 11.d0/4.d0-nf/6.d0;
        b1 = 51.d0/8.d0-19.d0*nf/24.d0;
        b2 = 2857.d0/128.d0 - nf*(5033.d0/1152.d0-nf*325d0/3456d0)
        b3 = 149753d0/1536d0+nf*nf*nf*1093d0/186624d0 +nf*nf*(50065d0+
     &       12944d0*Zeta3)/41472d0 -nf*(1078361d0 +39048d0*Zeta3)/
     &       41472d0 + 891d0*Zeta3/64d0

        nsteps = 100000

        steps = dlog(mur*mur/MZ/MZ)/(1.d0*nsteps)
c        write(*,*)'steps =',steps
  
        res = as_mz

C        RK4 method, without any threshold

        if (N.eq.0) then
        do i = 0 ,nsteps

        k1 = -steps*b0*res*res
        k2 = -steps*b0*(res+k1/2.0d0)*(res+k1/2.0d0)
        k3 = -steps*b0*(res+k2/2.0d0)*(res+k2/2.0d0)
        k4 = -steps*b0*(res+k3)*(res+k3)
        res = res + (k1+2d0*k2+2d0*k3+k4)/6.d0
        enddo
        elseif (N.eq.1) then
        do i = 0, nsteps
        k1=-steps*res*res*(b0 +res*b1)
        k2=-steps*(res+k1/2.d0)*(res+k1/2.d0)*(b0 +(res+k1/2.d0)*b1)
        k3=-steps*(res+k2/2.d0)*(res+k2/2.d0)*(b0 +(res+k2/2.d0)*b1)
        k4=-steps*(res+k3)*(res+k3)*(b0 +(res+k3)*b1)
        res = res + (k1+2*k2+2*k3+k4)/6.d0
        enddo 
        elseif (N.eq.2) then
        do i=0,nsteps
        k1 = -steps*res*res*(b0 +res*b1+res*res*b2)
        k2 = -steps*(res+k1/2.d0)*(res+k1/2.d0)*(b0 +
     &      (res+k1/2.d0)*b1 +(res+k1/2.d0)*(res+k1/2.d0)*b2);
        k3 = -steps*(res+k2/2.d0)*(res+k2/2.d0)*(b0 +(res+k2/2.d0)*b1 + 
     &      (res+k2/2.d0)*(res+k2/2.d0)*b2)
        k4 = -steps*(res+k3)*(res+k3)*(b0 +(res+k3)*b1 +
     &      (res+k3)*(res+k3)*b2)
        res = res + (k1+2*k2+2*k3+k4)/6.d0
        enddo
        elseif (N.eq.3) then
        do i=0,nsteps
        k1 = -steps*res*res*(b0 +res*b1 +res*res*b2 +res*res*res*b3)
        k2=-steps*(res+k1/2d0)*(res+k1/2d0)*(b0 +(res+steps*k1/2d0)*b1+
     &    (res+steps*k1/2.d0)*(res+k1/2.d0)*b2 +(res+steps*k1/2d0)*
     &    (res+k1/2d0)*(res+k1/2d0)*b3)
        k3 = -steps*(res+k2/2.d0)*(res+steps*k2/2.d0)*(b0 +
     &     (res+k2/2.d0)*b1 +(res+k2/2.d0)*(res+k2/2.d0)*b2 +
     &     (res+k2/2.d0)*(res+k2/2.d0)*(res+k2/2.d0)*b3)
        k4 = -steps*(res+k3)*(res+k3)*(b0 +(res+k3)*b1 +
     &     (res+k3)*(res+k3)*b2 + (res+k3)*(res+k3)*(res+k3)*b3)
        res = res + (k1+2d0*k2+2d0*k3+k4)/6.d0
        enddo
        endif

        as_n3loxs = res
c      write(*,*)'alphas =', mur,N,mz,as_mz,res

        return
        end

