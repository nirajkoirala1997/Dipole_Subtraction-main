      program uU2eE_Virtual 
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL
      character*50 name
      external flo2_Vir
      include 'coupl.inc'
      include 'nexternal.inc'
      call setpara('param_card.dat',.true.)
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

      print*,'  '
      print*,"Press 1 to initialise VEGAS:"
      print*,"Press 2 to initialise CUBA-VEGAS:"
        read*,i
        if (i .eq.  1) then
        Print*,"2. Calculating Virtual Contribution"
        print*,"````````````````````````````````````"
        Print*,"    Virtual + Dipole over 2-Body PS"
        print*," "
        print*," "

        call brm48i(40,0,0) ! initialize random number generator
        call vsup(3,npt1,its1,flo2_Vir,ai_lo2,sd,chi2)
c        ai_lo2=1.079582        
        print*,"  "
        write(*,*)'The answer is =', ai_lo2
        write(*,*)"Integral      =",ai_lo2,"+-",sd
        write(*,*)"with chisq    =",chi2
        print*,"  "
        print*,"  "
c
        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
