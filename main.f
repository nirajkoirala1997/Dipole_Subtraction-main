      program Drell_yan_dipoleSubtraction 
      implicit double precision (a-h,o-z)
      integer leg
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/leg_choice/leg
      common/usedalpha/AL
      common/set/set1
      common/countc/n4
      character*50 name
      external fnlo3
      external dipole_uU_g
      external flo2_PK
      external flo2_Vir
      include 'coupl.inc'
      include 'nexternal.inc'
      call setpara('param_card.dat',.true.)
        set1=0
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
      
c      am1 = 0.51099895000d-3
      am1=0.0d0
      am2=0.0d0
      am3=0d0
      am4=0d0
      am5=0d0
      leg=0
      ! energy
      s=ecm*ecm
      Al = 0.118d0

      n4 = 0

      print*,'  '
      print*,"Press 1 to initialise VEGAS:"
      print*,"Press 2 to initialise CUBA-VEGAS:"
        read*,i
        IF (I .EQ. 1) THEN
        print*,"----------------------------------"
        print*,"|Initializing Dipole Subtraction  |"
        print*,"----------------------------------"
        print*," "
        print*,"1. real - Dipole Over 3body phasespace"
        print*,"``````````````````````````````````````"
        print*," "
        print*," "
c        goto 1990
       call brm48i(40,0,0) ! initialize random number generator
       call vsup(6,npt1,its1,fnlo3,ai_nlo3,sd,chi2)
c       ai_nlo3=28.0150449
        print*,"  "
        write(*,*)'The answer is =', ai_nlo3
       write(*,*)"Integral      =",ai_nlo3,"+-",sd
       write(*,*)"with chisq    =",chi2
       write(*,*)"Unphysical count =",n4
        print*," "
        print*," "
        stop
1990        Print*,"2. Calculating Virtual Contribution"
        print*,"````````````````````````````````````"
        Print*,"    Virtual + Dipole over 2-Body PS"
        print*," "
        print*," "

c        call brm48i(40,0,0) ! initialize random number generator
c        call vsup(3,npt1,its1,flo2_Vir,ai_lo2,sd,chi2)
c        ai_lo2=1.9427896217405324E-003
        ai_lo2=1.079582        
        print*,"  "
        write(*,*)'The answer is =', ai_lo2
        write(*,*)"Integral      =",ai_lo2,"+-",sd
        write(*,*)"with chisq    =",chi2
        print*,"  "
        print*,"  "


        print*,"3. Calculating P and K terms for both Legs"
        print*,"````````````````````````````````````"
        print*," "
        print*," "
         PK=-0.409942883208528E-02
         print*,"Sum PK =",PK 
        print*," "
        print*," "
        print*," "
        print*," "
         print*,"Total Sigma :",ai_nlo3 + ai_lo2! + PK
c        leg=1
c        call brm48i(40,0,0) ! initialize random number generator
c        call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)
c        print*,"  "
c        write(*,*)'The answer is  =', ai_lo2
c        write(*,*)"Integral of PK1=",ai_lo2,"+-",sd
c        write(*,*)"with chisq    =",chi2
c        print*,"  "
c
c        print*,"4. Calculating P and K terms for Leg 2"
c        print*,"````````````````````````````````````"
c        leg=2
cc        call brm48i(40,0,0) ! initialize random number generator
cc        call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)
c        print*,"  "
c        write(*,*)'The answer is   =', ai_lo2
c        write(*,*)"Integral of PK2 =",ai_lo2,"+-",sd
c        write(*,*)"   with chisq   =",chi2
c        print*,"  "
c
        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
