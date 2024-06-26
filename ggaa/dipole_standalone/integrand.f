C Cuba Specific function.      
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
      integer function integrand(ndim,xx,ncomp,f,userdata,nvec ,core
     .                            ,weight,iter)
      implicit none 
c     cuba specific parameters
      integer n4,ndim,ncomp,nvec,core,iter,userdata 
      real*8 xx(ndim) ,f(ncomp),weight,fnlo3
      common/countc/n4
      external fnlo3
      f(1) = fnlo3(xx,weight)
      return
      end
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C        
c      Our vegas specific
      double precision function fnlo3(xx,weight)
      implicit double precision(a-h,o-z) 
C -------------------------------------------------------------------- C        
      dimension  xx(10),p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),q(0:3)
     .             ,p(0:3,1:5),dip(27),B1(1:2),xl(15),f1(-6:6),f2(-6:6)
     .             , sig(1:5),SumD(1:2)
      integer i35,i45,is5,itest
      integer k1,k2,k3,ipass,n4,unphy
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+6)
      common/amass/am1,am2,am3,am4,am5
      common/energy/s
      common/set/set1
      common/distribution/xq
      common/countc/n4
      common/usedalpha/AL,ge
      common/scales/xinvmass
      common/counter/ifilter,itot_ev,iselect_scale
      common/counter_diff/diff,eps
      common/t_cuts/e_cut,t_cut
c      common/momenta5/p1,p2,p3,p4,p5
      external dipole_gq_q




      xa = xx(1)
      xb = xx(2)
      rsp = dsqrt(xa*xb*s)

      xmz = 91.1876d0

      ipass1 = 0
       fnlo3 = 0

c      eps = 2.0d0
      xlow = xq - eps
      xhigh = xq + eps




        call kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)
        call cuts3(p1,p2,p3,p4,p5,rsp,icuts,inf_PS)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Technical Cut

       soft  = e_cut
       coll1 = t_cut
       coll2 = t_cut 
 
       i15=0d0
       i25=0d0
       is5=0d0
       itest=0d0

        s15 = 2d0*dot(p1,p5)
        s25 = 2d0*dot(p2,p5)
        s12=am1**2 + am2**2 + 2d0*dot(p1,p2)
        s35=am3**2 + am5**2 + 2d0*dot(p3,p5)
        s45=am4**2 + am5**2 + 2d0*dot(p4,p5)
        s34=am3**2 + am4**2 + 2d0*dot(p3,p4)



c     collinear
       if (dabs(s15) .lt. coll1) i15=1
       if (dabs(s25) .lt. coll2) i25=1

c      soft
       e5=0.5d0*(s12-s34)/rsp
 
       if(e5 .le. soft) is5=1
        itest=is5+i15+i25
c        unphy = unphy + itest



c        print*,s15,s25,e5
c        print*,coll1,coll2,soft
c        print*," "
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c        if (unphy .eq. 0 .and. icuts .eq. 1 .and. itest=0) then ! with zero unphysical PS points proceed
        if (unphy .eq. 0 .and. icuts .eq. 1 ) then ! with zero unphysical PS points proceed
     
c        itot_ev = itot_ev + 1

        scale = xinvmass
        if ( scale .ge. xlow .and. scale .le. xhigh) then
c         if (inf_PS .eq. 1 ) goto 151

c         iselect_scale = iselect_scale + 1 

          xmuf=scale
          xmur=scale

          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call setlum(f1,f2,xl)

          AL = alphasPDF(xmur)

          call p1dtop2d_5(p1,p2,p3,p4,p5,p)
          call amp_mat_r(p,sig)


          SumD(1) = dipole_gg_g(1,p) + dipole_gg_g(2,p)

          sigma = xl(4)*( sig(4)-SumD(1) )

c          if (sig(4) - SumD(1) .ge. diff) then
c                  ifilter = ifilter + 1 
c           goto 151
c          endif
c          if (sig(4) .ge. 100d0) print*,"|M^2|:",sig(4),"SumD:",SumD(1) 
c
 


c          sigma = xl(4)*sig(4)
c          sigma = xl(4)*SumD(1)

c          call dipole(p,dip,sum_dipole)
c          print*,"sum","sigma",sum_dipole,sigma
c          stop
c          sigma = xl(4)*SumD(1)
c          sigma = xl(4)*sum_dipole


          pi_1 = 0.5d0*rsp
          flux = 4d0*pi_1*rsp
          xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
          wgt=xxjac*xnorm*sigma*weight
          fnlo3=wgt/weight/2d0/eps
          endif
c        print*,sigma,sig(4),dipole_gg_g(1,p),dipole_gg_g(2,p),weight
          endif

c        if (fnlo3 .ne. fnlo3) then
c        print*,sigma,sig(4),dipole_gg_g(1,p),dipole_gg_g(2,p),weight
cc        stop
c        endif

151   return
      end
c---------------------------------------------------------------------
