      program uU2eE_Virtual 
      implicit double precision (a-h,o-z)
      dimension x(10),y(10),ai_lo2(1:50),err(0:50)
      parameter (pi=3.14159265358979d0)

      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL,ge   
      common/distribution/xq
      character*50 name
      character*100 run_tag,filename
      external flo2_LO

      !input data card
      open(unit=10,file='run.vegas.dat',status='unknown')    
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
      close(10)

      open(unit=10,file='param_card.dat',status='unknown')    
      read (10,*) ge       ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          !lhapdf set
      read (15,*) it_max        !lhapdf set
      read (15,*) xq_initial
      read (15,*) step_size         !step
      read (15,*) run_tag           !run_tag saves the output 
      read (15,*) iprint            !save data in output file ../summary
      close(15)


        call initpdfsetbyname(name)
        Call initPDF(0)
      
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
        if (i .eq. 1 ) then
        print*," "
        print*," "
        print*,"____________________________________"
        Print*," Calculating LO_DY"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        print*," "
        print*," "
        xq = xq_initial

c       writes data in output file
        filename = "LO.dat"
        if (iprint .eq. 1) call output(run_tag,filename)
        
        do j=1,it_max

          print*," "
      write(*,*) achar(27)//'[1;33m' // "For xq = ",int(xq) ,achar(27)
     .   //'[0m'
          print*," "

         call brm48i(40,0,0) ! initialize random number generator
         call vsup(3,npt1,its1,flo2_LO,ans,sd,chi2)
         ai_lo2(j) = ans
         err(j) = sd

            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral=",
     .  ai_lo2(j),achar(27) //'[0m', "+-",err(j)
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "
            xq = xq + step_size
        enddo

           xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq               LO Integral",
     .   "                   error",
     . achar(27)//'[0m'
        do j=1,it_max
          write(*,'(i7,3e27.15)')int(xq),ai_lo2(j),err(j)
          xq = xq + step_size
        enddo

        if (iprint .eq. 0) goto 123
        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .  //trim(filename),status='unknown')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_lo2(i),err(i)
          xq = xq + step_size
         enddo
         close(20)
123         continue
c
        elseif(I .eq. 2) THEN
c                CALL cubacheck
        endif
       end
