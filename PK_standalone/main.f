      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name
      common/energy/s
      external flo2_PK
      common/leg_choice/leg
      common/usedalpha/AL
      common/distribution/xq
      

      include 'coupl.inc'
      include 'nexternal.inc'
      call setpara('param_card.dat',.true.)

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          !lhapdf set
      read (15,*) iorder        !iorder no of q for distribution
      read (15,*) xq            ! initialise xq value
      close(15)
        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(1)
       s=ecm*ecm
       print*,"  "
       print*,"  "
        print*,"____________________________________"
        print*,"1. Calculating P and K terms for Leg 1"
        print*,"____________________________________"
        print*,"````````````````````````````````````"

        leg=1
        do j=0,iorder
          print*,"For xq:",xq_initial
          call brm48i(40,0,0) ! initialize random number generator
          call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)

          xq=xq+50d0
          PK1 = ai_lo2

          print*,"  "
          write(*,*)'The answer is  =', ai_lo2
          write(*,*)"Integral of PK1=",ai_lo2,"+-",sd
          write(*,*)"with chisq    =",chi2
          print*,"  "
        enddo

        print*,"____________________________________"
        print*,"2. Calculating P and K terms for Leg 2"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        leg=2
        xq = xq_initial
        print*," "
        print*,"For xq:",xq
        do j=0,iorder
c          call brm48i(40,0,0) ! initialize random number generator
c          call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)
          PK2= ai_lo2
          xq=xq+50d0
          print*,"  "
          write(*,*)'The answer is   =', ai_lo2
          write(*,*)"Integral of PK2 =",ai_lo2,"+-",sd
          write(*,*)"   with chisq   =",chi2
          print*,"  "
        enddo  
c          write(*,'(A,3e27.15)')"Total PK1 + PK2 =",PK1+PK2
cc          print*,"  "
c          print*,"  "


       end
      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end
