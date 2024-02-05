
      program Drell_yan_dipoleSubtraction 
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      character*50 name
      external fnlo3
      external dipole_uU_g
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
      ! energy
      s=ecm*ecm
      print*,'  '
      print*,"Press 1 to initialise VEGAS:"
      print*,"Press 2 to initialise CUBA-VEGAS:"
        read*,i
        IF (I .EQ. 1) THEN
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(6,npt1,its1,fnlo3,ai_nlo3,sd,chi2)
        print*,"  "
        write(*,*)'The answer is =', ai_nlo3
        write(*,*)"Integral      =",ai_nlo3,"+-",sd
        write(*,*)"with chisq    =",chi2
        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
