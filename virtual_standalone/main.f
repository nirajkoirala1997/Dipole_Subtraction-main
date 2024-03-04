      program uU2eE_Virtual 
      implicit double precision (a-h,o-z)
      dimension x(10),y(10),ai_lo2(1:50)
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL,Al_ew
      common/distribution/xq
      character*50 name
      character*100 run_tag,filename
      external flo2_Vir

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')    
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1*2
      close(10)


      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) Al_ew       ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          !lhapdf set
      read (15,*) it_max        !lhapdf set
      read (15,*) xq_initial
      read (15,*) xstep         !step
      read (15,*) run_tag       ! dir name to save data
      read (15,*) iprint            !save data in output file ../summary
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
c
c        
c        open(unit=17,file='Output.dat',status='unknown')
c        do i=1,6
c        read(17,*) y(i)
c        enddo
c        close(17)
c        do i=1,6
c        print*,x(i)/y(i)
c        enddo
c        stop

        i=1
        if (i .eq.  1) then
        print*," "
        print*," "
        print*,"____________________________________"
        Print*,"    Calculating Virtual Contribution"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        Print*,"    Virtual + Dipole over 2-Body PS"
        print*," "
        print*," "
        xq = xq_initial
c       writes data in output file
        filename = "virtual.dat"
        if(iprint .eq. 1)  call output(run_tag,filename)

        do j=1,it_max

          print*," "
      write(*,*) achar(27)//'[1;33m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "
c        print*,"To write result in Outputfile press 1 else 2:"
c        read*,j
         call brm48i(40,0,0) ! initialize random number generator
         call vsup(3,npt1,its1,flo2_Vir,ai_lo2(j),sd,chi2)
c         if (j .eq. 1) then
c         open(unit=17,file='../Output_low_q_LO.dat',status='unknown',
c     .          position='append')
c         write(17,*) xq,ai_lo2
c         close(17)
c         endif
         
         xq = xq + xstep

            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral=", 
     .  ai_lo2(j),achar(27) //'[0m', "+-",sd
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "

        enddo

        xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
     . achar(27)//'[0m'
        do j=1,it_max
          write(*,'(i7,3e27.15)')
     .             int(xq),ai_lo2(j)
          xq = xq + xstep
        enddo
        if(iprint .eq. 0) goto 123
        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .  //trim(filename),status='unknown')
c     .  //trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_lo2(i)
          xq = xq + xstep
         enddo
         close(20)
123     continue 
c
        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
