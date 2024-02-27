      program uU2eE_Virtual 
      implicit double precision (a-h,o-z)
      dimension x(10),y(10),ai_lo2(1:50)
      parameter (pi=3.14159265358979d0)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL
      common/distribution/xq
      character*50 name
      character*100 run_tag
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
      read (15,*) it_max        !lhapdf set
      read (15,*) xq_initial
      read (15,*) xstep         !step
      read (15,*) run_tag       ! dir name to save data
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
        Print*,"3. Calculating Virtual Contribution"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        Print*,"    Virtual + Dipole over 2-Body PS"
        print*," "
        print*," "
        xq = xq_initial
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

           call output(run_tag)
        open(unit=20,file='../summary/'//trim(run_tag)//'/virtual.dat',
     .                     status='unknown')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_lo2(i)
          xq = xq + xincr
         enddo
         close(20)
       
c
        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
       
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
        subroutine output(run_tag)
        implicit none 
        integer ierr
        character*100 dir_path,run_tag

       ! Check if the directory exists
       dir_path ="../summary/"// trim(run_tag)   
       call system("test -d "// dir_path // " && echo 1 > 
     .   command.txt || echo 0  >
     .        command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr
       close(13)
         print*," "
        if (ierr .eq. 1) print*,"Directory exists overwriting data in"
       call system("rm command.txt")
       if( ierr .ne. 1) then 
         print*,"Directory not found " //dir_path 
         print*,"Making new directory.."
         print*," "
         print*,"Writing Data in"
         CALL SYSTEM("cd ../summary && mkdir -p " // dir_path) 
       endif
         CALL SYSTEM("cd " //dir_path// 
     . " && echo $PWD")
        print*," "
        end
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
