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
      read (15,*) iprint        ! to print data in file
      close(15)


        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(1)
       s=ecm*ecm
       print*,"Enter Leg Choice 1 or 2"
       read*,leg_user




c~~~~~~~~~~~~~~~~~~~~~~~~~~ Leg 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~c       

       if (leg_user .eq. 1 ) then
        print*,"  "
        print*,"  "
        print*,"____________________________________"
        print*,"Calculating P and K terms for Leg 1"
        print*,"____________________________________"
        print*,"````````````````````````````````````"

        leg=1
        filename = 'PK1.dat'
c ~~~~~~~~~~~Writing in a file to store~~~~~~~~~~~~c        
      if (iprint .eq. 1) call output(run_tag,filename)            

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
c~~~~~~~~~~~~~~~~~~~~~~~~~ * END * ~~~~~~~~~~~~~~~~~~~~~~~~~~c       





c~~~~~~~~~~~~~~~~~~~~~~~~~~ Leg 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~c       
      elseif (leg_user .eq. 2) then
        print*,"____________________________________"
        print*,"Calculating P and K terms for Leg 2"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        print*," "
        leg=2
      filename = 'PK2.dat'
c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
      if (iprint .eq. 1) call output(run_tag,filename)            
c ~~~~~~~~~~~~uncomment this to write~~~~~~~~~~~~~~~~c
        xq = xq_initial

        do j=1,it_max
          print*," "
      write(*,*) achar(27)//'[1;33m' // "For xq = ",int(xq) ,achar(27) 
     .   //'[0m'
          print*," "
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call brm48i(40,0,0) ! initialize random number generator 
          call vsup(4,npt1,its1,flo2_PK,ai_lo2,sd,chi2)            
c___________________________________________________________________         

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




c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .eq. 0 ) goto 123
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

        if (iprint .eq. 0 ) goto 123
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
123        continue
       end
      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end
