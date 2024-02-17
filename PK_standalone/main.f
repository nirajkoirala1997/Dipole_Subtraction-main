      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name
      common/energy/s
      external flo2_PK
      common/leg_choice/leg

      include 'coupl.inc'
      include 'nexternal.inc'
      call setpara('param_card.dat',.true.)

      name='CT10nlo'
      call initpdfsetbyname(name)
      Call initPDF(1)

      !input data card
      open(unit=10,file='run.vegas.dat',status='unknown')
      read (10,*) npt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      close(10)

      open(unit=15,file='run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      close(15)
        s=ecm*ecm
       print*,"  "
       print*,"  "
        print*,"____________________________________"
        print*,"1. Calculating P and K terms for Leg 1"
        print*,"____________________________________"
        print*,"````````````````````````````````````"

        leg=1
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)
        PK1 = ai_lo2
       print*,"  "
        write(*,*)'The answer is  =', ai_lo2
        write(*,*)"Integral of PK1=",ai_lo2,"+-",sd
        write(*,*)"with chisq    =",chi2
        print*,"  "

        print*,"____________________________________"
        print*,"2. Calculating P and K terms for Leg 2"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        leg=2
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)
        PK2= ai_lo2
        print*,"  "
        write(*,*)'The answer is   =', ai_lo2
        write(*,*)"Integral of PK2 =",ai_lo2,"+-",sd
        write(*,*)"   with chisq   =",chi2
        print*,"  "
        write(*,'(A,3e27.15)'),"Total PK1 + PK2 =",PK1+PK2
        print*,"  "
        print*,"  "


       end
      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end
