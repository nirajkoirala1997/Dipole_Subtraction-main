      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_PK,flo2_PKDel,flo2_PKReg
      common/leg_choice/leg
      common/usedalpha/AL,ge   
      common/distribution/xq
      dimension PKPlus(1:50),err_Plus(1:50)
      dimension PKReg(1:50),err_Reg(1:50)
      dimension PKDel(1:50),err_Del(1:50)
      dimension PK(1:50),err(1:50)
      

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')
      do i=1,6
      read (10,*)
      enddo
      read (10,*) pt1           ! vegas points     
      read (10,*) its1          ! vegas iterations 
      npt1 = pt1
      close(10)
      

      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) ge          ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          ! lhapdf set
      read (15,*) it_max        ! it_max no of q for distribution
      read (15,*) xq            ! initialise xq value
      read (15,*) xincr         ! increment in Gev from xq 
      read (15,*) run_tag
      read (15,*) iprint        ! to print data in file
      close(15)

c ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        

      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*)
      read (20,*) 
      read (20,*) filename
      close(20)
      if (iprint .eq. 1) call output(run_tag,filename)            
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        



        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(0)
       s=ecm*ecm

        print*,"  "
        print*,"____________________________________"
        print*,"    Calculating P and K terms       "
        print*,"____________________________________"
        print*,"````````````````````````````````````"



c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Plus  functions ]
        print*,"  "
        print*,"____________________________________"
        print*,"  Calculating Regular and [+] terms "
        print*,"____________________________________"

        do j=1,it_max

          print*," "
      write(*,*) achar(27)//'[1;32m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

          call brm48i(40,0,0) ! initialize random number generator
          call vsup_2(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)

          PKPlus(j)   = ai_lo2
          err_plus(j) =sd
            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral Plus ", 
     .      PKPlus(j),achar(27) //'[0m', "+-",err_plus(j)
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "

           xq=xq + xincr
        enddo
        xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","   Integral Plus",
     .  "                    error",
     . achar(27)//'[0m'
        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKPlus(j),err_plus(j)
          xq = xq + xincr
        enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ regular functions ]
        xq = xq_initial
        print*,"  "
        print*,"____________________________________"
        print*,"    Calculating Regular terms"
        print*,"____________________________________"

        do j=1,it_max

          print*," "
      write(*,*) achar(27)//'[1;32m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

          call brm48i(40,0,0) ! initialize random number generator
          call vsup(4,npt1,its1,flo2_PKReg,ai_lo2,sd,chi2)

          PKReg(j) = ai_lo2
          err_Reg(j)=sd
            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral Regular", 
     .      ai_lo2,achar(27) //'[0m', "+-",err(j)
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "

           xq=xq + xincr
        enddo
        xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","   Integral Regular",
     .  "                    error",
     . achar(27)//'[0m'
        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKReg(j),err_Reg(j)
          xq = xq + xincr
        enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ delta functions ]
        xq = xq_initial
        print*,"____________________________________"
        print*," Calculating Delta function terms"
        print*,"____________________________________"

        do j=1,it_max

          print*," "
      write(*,*) achar(27)//'[1;32m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "

          call brm48i(40,0,0) ! initialize random number generator
          call vsup(3,npt1,its1,flo2_PKDel,ai_lo2,sd,chi2)

          PKDel(j) = ai_lo2
          err_Del(j)=sd
            print*,"  "
            print*,"  "
            write(*,*)achar(27)//'[1;32m'//"Integral Delta ", 
     .      ai_lo2,achar(27) //'[0m', "+-",sd
            write(*,*)"with chisq    =",chi2
            print*," "
            print*," "

           xq=xq + xincr
        enddo
        xq = xq_initial
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","  Integral PK ",
     .  "                    error",
     . achar(27)//'[0m'

        do j=1,it_max
        PK(j) = PKPlus(j) + PKReg(j) + PKDel(j)
        err(j) = err_Plus(j) + err_Reg(j) + err_Del(j)
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PK(j),err(j)
          xq = xq + xincr
        enddo
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[  * END * ]      

c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .eq. 0 ) goto 123
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
c     .   '/'//trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(21,*)xq,PK(i),err(i)
          xq = xq + xincr
         enddo
         close(21)
123        continue
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        


      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end
