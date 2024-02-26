      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name
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
      close(15)
        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(1)
       s=ecm*ecm
       print*,"Enter Leg Choice 1 or 2"
       read*,leg_user
       if (leg_user .eq. 1 ) then
        print*,"  "
        print*,"  "
        print*,"____________________________________"
        print*,"Calculating P and K terms for Leg 1"
        print*,"____________________________________"
        print*,"````````````````````````````````````"

        leg=1
         
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
      elseif (leg_user .eq. 2) then
        print*,"____________________________________"
        print*,"Calculating P and K terms for Leg 2"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        print*," "
        leg=2
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
      endif
        xq = xq_initial
        print*,"Total PK1 + PK2 for" 

       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
     . achar(27)//'[0m'
       print*," "
        do j=1,it_max
          xq = xq_initial*j
          write(*,'(i7,3e27.15)')
     .             int(xq),PK1(j)+PK2(j)
        enddo


       end
      double precision function dot(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end
