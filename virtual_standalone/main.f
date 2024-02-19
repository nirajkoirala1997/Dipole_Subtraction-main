      program uU2eE_Virtual 
      implicit double precision (a-h,o-z)
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL
      common/distribution/xq
      character*50 name
      external flo2_Vir
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
      read (15,*) iorder        !lhapdf set
      close(15)

        call initpdfsetbyname(name)
        Call initPDF(1)
      
c      am1 = 0.51099895000d-3
      am1=0.0d0
      am2=0.0d0
      am3=0d0
      am4=0d0
      am5=0d0
      leg=0
      ! energy
      s=ecm*ecm

c      print*,'  '
c      print*,"Press 1 to initialise VEGAS:"
c      print*,"Press 2 to initialise CUBA-VEGAS:"
c        read*,i
        i=1
        if (i .eq.  1) then
        print*," "
        print*," "
        print*,"____________________________________"
        Print*,"3. Calculating Virtual Contribution"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        Print*,"    Virtual + Dipole over 2-Body PS"
        print*," "
        print*," "
        xq=100d0
        do i=0,iorder

        print*," "
        print*,"For xq:",xq
        print*," "

         call brm48i(40,0,0) ! initialize random number generator
         call vsup(3,npt1,its1,flo2_Vir,ai_lo2,sd,chi2)
         
         xq = xq + 50d0
c        ai_lo2=1.079582        
         print*,"  "
         write(*,*)'The answer is =', ai_lo2
         write(*,*)"Integral      =",ai_lo2,"+-",sd
         write(*,*)"with chisq    =",chi2
         print*,"  "
         print*,"  "

        enddo
c
        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
