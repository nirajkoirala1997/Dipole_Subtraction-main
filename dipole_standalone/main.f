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
        name='CT10nlo'
        call initpdfsetbyname(name)
        Call initPDF(1)

      !input data card
      open(unit=10,file='run.vegas.dat',status='unknown')    
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
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
c        i=1
        IF (I .EQ. 1) THEN
        print*,"----------------------------------"
        print*,"|Initializing Dipole Subtraction  |"
        print*,"----------------------------------"
        print*," "
        print*,"1. real - Dipole Over 3body phasespace"
        print*,"``````````````````````````````````````"
        print*," "
        print*," "
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

        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
