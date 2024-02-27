      program Drell_yan_dipoleSubtraction 
      implicit double precision (a-h,o-z)
      integer leg
      character*50 green
      parameter (pi=3.14159265358979d0)
      dimension ai_nlo3(0:50)
      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/leg_choice/leg
      common/usedalpha/AL
      common/set/set1
      common/countc/n4
      common/distribution/xq
      character*50 name
      character*100 run_tag,filename
      external fnlo3
      external dipole_uU_g
      green = ''//achar(27)//'[32m'//char(27)//'[0m'

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')    
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable 
      read (15,*) run_tag               ! name of run directory to save output
      close(15)
        name = 'MMHT2014nlo68cl' 
      
        call initpdfsetbyname(name)
        Call initPDF(1)
c        Al = alslhaPDF(mur)
c      am1 = 0.51099895000d-3
      am1=0.0d0
      am2=0.0d0
      am3=0d0
      am4=0d0
      am5=0d0
      leg=0
      ! energy
      s=ecm*ecm
      filename = "real.dat"
      print*,'  '
c      print*,"Press 1 to initialise VEGAS:"
c      print*,"Press 2 to initialise CUBA-VEGAS:"
c        read*,i
        i=1
        IF (I .EQ. 1) THEN
          print*,"----------------------------------"
          print*,"|Initializing Dipole Subtraction  |"
          print*,"----------------------------------"
          print*," "
          print*,"4. real - Dipole Over 3body phasespace"
          print*,"``````````````````````````````````````"
          print*," "
          print*," "

c  saves the data in output file           
          call output(run_tag,filename)

          xq = xq_initial
          do j = 1,it_max
            
c            print*,"For xq=",xq ! base value
c            write(*, '(A)') // green // "For xq = " // reset //
          print*," "
      write(*,*) achar(27)//'[1;33m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

            call brm48i(40,0,0) ! initialize random number generator
            call vsup(6,npt1,its1,fnlo3,ai_nlo3(j),sd,chi2)

            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral=", 
     .  ai_nlo3(j),achar(27) //'[0m', "+-",sd
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "
            xq = xq + step_size 
          enddo
          
          xq = xq_initial
       write(*,*)achar(27)//'[1;92m'//"   xq"," ","           Integral",
     . achar(27)//'[0m'
          
          do j=1,it_max
          write(*,'(i7,3e27.15)')int(xq),ai_nlo3(j)
          
          xq = xq + step_size 
          enddo

        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .          //trim(filename),
     .                     status='unknown')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_nlo3(i)
          xq = xq + step_size 
         enddo
         close(20)
        elseif(I .eq. 2) THEN
                CALL cubacheck
        endif
       end
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        

c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
        subroutine output(run_tag,filename1)
        implicit none 
        integer ierr1,ierr2
        character*100 dir_path,run_tag,dir_pathtmp,decision
     . ,filename1,filename

       ! Check if the directory exists
       dir_path ="../summary/"// trim(run_tag)
       filename ="../summary/"// trim(run_tag) //"/"// trim(filename1)  
       dir_pathtmp ="../summary/temp_"// trim(run_tag)   
c checking Directory
        call system("test -d "// dir_path // " && echo 1 > 
     .   command.txt || echo 0  >
     .        command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr1
       close(13)
       call system("rm command.txt")
       if( ierr1.ne. 1) then 
         print*,"Directory not found " //dir_path 
         print*,"Making new directory.."
         print*," "
         CALL SYSTEM("cd ../summary && mkdir -p " // dir_path) 
       endif
c checking file
       call system("test -f "// filename // " && echo 1 > 
     .   command.txt || echo 0  >
     .        command.txt")
       open(unit=13,file="command.txt",status="unknown")
       read(13,*)ierr2
       close(13)
       call system("rm command.txt")
c Proceed with decision       

       if (ierr2 .eq. 1) then
        print*,"File  already exist  "//filename
        print*,"All previous data for this run will be lost. Overwrite ?
     .     [y/n]"
        read*,decision
        if (decision .eq. 'y') then
        print*,"Overwriting data in"// filename
        else
           stop
        endif
        endif
        call sleep(3)
        end
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
