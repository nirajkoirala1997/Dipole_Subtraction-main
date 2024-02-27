      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_PK
      common/leg_choice/leg
      common/usedalpha/AL
      common/distribution/xq
      dimension PK1(1:50),PK2(1:50)
      

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
      read (15,*) name          ! lhapdf set
      read (15,*) it_max        ! it_max no of q for distribution
      read (15,*) xq            ! initialise xq value
      read (15,*) xincr         ! increment in Gev from xq 
      read (15,*) run_tag
      close(15)


        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(1)
       s=ecm*ecm
       print*,"Enter Leg Choice 1 or 2"
       read*,leg_user
c~~~~~~~~~~~~~~~~ Leg 1~~~~~~~~~~~~~~~~c       
       if (leg_user .eq. 1 ) then
        print*,"  "
        print*,"  "
        print*,"____________________________________"
        print*,"Calculating P and K terms for Leg 1"
        print*,"____________________________________"
        print*,"````````````````````````````````````"

        leg=1
      filename = 'PK1.dat'
c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        call output(run_tag,filename)            
c ~~~~~~~~~~~~proceed to write~~~~~~~~~~~~~~~~c
       
         
        do j=1,it_max
          print*," "
      write(*,*) achar(27)//'[1;33m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

          call brm48i(40,0,0) ! initialize random number generator
          call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)

          PK1(j) = ai_lo2
            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral=", 
     .      ai_lo2,achar(27) //'[0m', "+-",sd
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "

           xq=xq + xincr
        enddo
        xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","       Integral PK1",
     . achar(27)//'[0m'
        do j=1,it_max
          write(*,'(i7,3e27.15)')
     .             int(xq),PK1(j)
          xq = xq + xincr
        enddo


c~~~~~~~~~~~~~~~ * END * ~~~~~~~~~~~~~~~~c       

c~~~~~~~~~~~~~~~~ Leg  2 ~~~~~~~~~~~~~~~~c       
      elseif (leg_user .eq. 2) then
        print*,"____________________________________"
        print*,"Calculating P and K terms for Leg 2"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        print*," "
        leg=2
      filename = 'PK2.dat'
c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        call output(run_tag,filename)            
c ~~~~~~~~~~~~proceed to write~~~~~~~~~~~~~~~~c
        xq = xq_initial

        do j=1,it_max
          print*," "
      write(*,*) achar(27)//'[1;33m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

          call brm48i(40,0,0) ! initialize random number generator
          call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)

          PK2(j) = ai_lo2

            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral=", 
     .  ai_lo2,achar(27) //'[0m', "+-",sd
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "
          xq=xq + xincr
        enddo  
        xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","       Integral PK2",
     . achar(27)//'[0m'
        do j=1,it_max
          write(*,'(i7,3e27.15)')
     .             int(xq),PK2(j)
          xq = xq + xincr
        enddo
      endif
c~~~~~~~~~~~~~~~ * END * ~~~~~~~~~~~~~~~~c       

c        xq = xq_initial
c        print*,"Total PK1 + PK2 for" 
c
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
c     . achar(27)//'[0m'
c       print*," "
c        do j=1,it_max
c          xq = xq + xincr
c          write(*,'(i7,3e27.15)')
c     .             int(xq),PK1(j)+PK2(j)
c        enddo


        if (leg_user .eq. 1) then
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
         xq = xq_initial
         do i=1,it_max
          write(21,*)xq,PK1(i)
          xq = xq + xincr
         enddo
         close(21)
        endif
c ~~~~~~~~~~~           ***              ~~~~~~~~~~~~c        
c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        

        if (leg_user .eq. 2) then
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
         xq = xq_initial
         do i=1,it_max
          write(21,*)xq,PK2(i)
          xq = xq + xincr
         enddo
         close(21)
        endif
c ~~~~~~~~~~~           ***              ~~~~~~~~~~~~c        


       end
      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end

cc ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
cc ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
c        subroutine output(run_tag)
c        implicit none 
c        integer ierr
c        character*100 dir_path,run_tag
c
c       ! Check if the directory exists
c       dir_path ="../summary/"// trim(run_tag)   
c       call system("test -d "// dir_path // " && echo 1 > 
c     .   command.txt || echo 0  >
c     .        command.txt")
c       open(unit=13,file="command.txt",status="unknown")
c       read(13,*)ierr
c       close(13)
c         print*," "
c        if (ierr .eq. 1) print*,"Directory exists overwriting data in"
c       call system("rm command.txt")
c       if( ierr .ne. 1) then 
c         print*,"Directory not found " //dir_path 
c         print*,"Making new directory.."
c         print*," "
c         print*,"Writing Data in"
c         CALL SYSTEM("cd ../summary && mkdir -p " // dir_path) 
c       endif
c         CALL SYSTEM("cd " //dir_path// 
c     . " && echo $PWD")
c        print*," "
c
c        end
cc ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
cc ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
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
        call sleep(2)
        end
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
c ~~~~~~~~~~~  *********************     ~~~~~~~~~~~~c        
