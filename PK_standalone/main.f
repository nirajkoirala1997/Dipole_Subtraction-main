      program intPK
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name,mode,mode1
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_Plus,flo2_PKDel,flo2_PKReg
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


c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[P and K terms from Here ] 
        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(0)
       s=ecm*ecm

      mode = "P and K terms"
      call printframe0(mode)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Plus  functions ]
      mode1 = "[+] distribution"
      pt2 = pt1*10
      npt2 = pt2
      its2 = its1

      call printframe0(mode1)
      call printframe1(pt2,its2)   ! Prints Vegas points

        do j=1,it_max

      call printframe2(xq)

c      -------------------------------------------------
         call brm48i(40,0,0) 
         call vsup(4,npt2,its2,flo2_Plus,ai_lo2,sd,chi2)
c      -------------------------------------------------

         PKPlus(j)   = ai_lo2
         err_plus(j) = sd

      mode = "Plus"
      call printframe3(mode,ai_lo2,sd,chi2)   

        xq=xq + xincr
        enddo

        xq = xq_initial
      call printframe4(mode)

        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKPlus(j),err_plus(j)
          xq = xq + xincr
        enddo
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ regular functions ]
        xq = xq_initial

      mode = "Regular Terms "
      call printframe0(mode)
      call printframe1(pt1,its1)

        do j=1,it_max

      call printframe2(xq)
c     -------------------------------------------------
      call brm48i(40,0,0) 
      call vsup(4,npt1,its1,flo2_PKReg,ai_lo2,sd,chi2)
c     -------------------------------------------------

          PKReg(j) = ai_lo2
          err_Reg(j)=sd

      mode = "Regular"
         call printframe3(mode,ai_lo2,sd,chi2)

           xq=xq + xincr
        enddo
        xq = xq_initial
      call printframe4(mode)

        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKReg(j),err_Reg(j)
          xq = xq + xincr
        enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ delta functions ]
        xq = xq_initial

      mode = "Delta Functions"
      call printframe0(mode)
      call printframe1(pt1,its1)

        do j=1,it_max

      call printframe2(xq)

c     -------------------------------------------------
        call brm48i(40,0,0) 
        call vsup(3,npt1,its1,flo2_PKDel,ai_lo2,sd,chi2)
c     -------------------------------------------------

          PKDel(j) = ai_lo2
          err_Del(j)=sd

      mode = "Delta"
      call printframe3(mode,ai_lo2,sd,chi2)

           xq=xq + xincr
        enddo
        xq = xq_initial

       call printframe4(mode)
        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKDel(j),err_Del(j)
          xq = xq + xincr
        enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Combining All ]
        xq = xq_initial
      print*,"  "
      write(*,*)achar(27)//'[1;32m'//
     . "   xq              Plus                      regular     
     .            delta                    combined PK     
     .    error",achar(27) //'[0m'

        do j=1,it_max
        PK(j) = PKPlus(j) + PKReg(j) + PKDel(j)
        err(j) = err_Plus(j) + err_Reg(j) + err_Del(j)
          write(*,'(i7,3e27.15,3e27.15,3e27.15,3e27.15,3e27.15)')
     .    int(xq),PKPlus(j),PKReg(j),PKDel(j),PK(j),err(j)
          xq = xq + xincr
        enddo
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[  * END * ]      


c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .ne. 1 ) goto 123
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


