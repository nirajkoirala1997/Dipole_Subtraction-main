      program intPK_Plus
      implicit double precision (a-h,o-z)
      dimension c(1:2)
      character*50 name,mode,mode1
      character*100 command,run_tag,dir_path,filename,filename1
      common/energy/s
      external flo2_Plus,flo2_PKDel,flo2_PKReg
      common/usedalpha/AL,ge   
      common/distribution/xq
      common/plus_cutoff/delta
      dimension PKPlus(1:50),err_Plus(1:50)
      dimension PKReg(1:50),err_Reg(1:50)
      dimension PKDel(1:50),err_Del(1:50)
      dimension PK(1:50),err(1:50)
      

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')
      do i=1,12     
        read (10,*)
      enddo
      read (10,*) pt1           ! vegas points     
      read (10,*) its1          ! vegas iterations 
      npt1 = pt1
      nptd = pt1/1000d0
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
      do i = 1,5
      read (20,*)
      enddo
      read (20,*) filename
      close(20)
      if (iprint .eq. 1) call output(run_tag,filename)            
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        

        delta  = 1d-10
        xq_initial = xq
      call initpdfsetbyname(name)
      Call initPDF(0)
       s=ecm*ecm

c      mode = "P and K terms"
c      call printframe0(mode)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ Plus  functions ]
      mode1 = "😇[+] distribution new implementation😇"

      call printframe0(mode1)
      print*,"Delta_cutoff=",delta
      call printframe1(pt1,its1)   ! Prints Vegas points

        do j=1,it_max

         call printframe2(xq)

c      -------------------------------------------------
         call brm48i(40,0,0) 
         call vsup(4,npt1,its1,flo2_Plus,ai_lo2,sd,chi2)
c      -------------------------------------------------


c     -------------------------------------------------
         call brm48i(40,0,0) 
         call vsup(2,nptd,its1,flo2_PKDel,ai_lo2_d,sd_d,chi2)
c     -------------------------------------------------

c        [Note-]
c        delta part contains finite pieces from 0 to 1-delta however
c        the singularity @x=1 is analytically removed with the 
c        singularity @x=1 from the first part as we are integrating 
c        from 0 to 1-delts hence delta part 
c        should have negative contribution to the integral.

         PKPlus(j)   = ai_lo2 - ai_lo2_d
         err_plus(j) = sd + sd_d !error cannot be subtracted

         mode = "Complete Plus"
         call printframe3(mode,PKPlus(j),err_plus(j),chi2)

         xq=xq + xincr
        enddo

        xq = xq_initial
        call printframe4(mode)

        do j=1,it_max
          write(*,'(i7,3e27.15,3e27.15)')
     .             int(xq),PKPlus(j),err_plus(j)
          xq = xq + xincr
        enddo

        xq = xq_initial

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[  * END * ]      


c ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
        if (iprint .ne. 1 ) goto 123
       open(unit=21,file='../summary/'//trim(run_tag)//
     .   '/'//trim(filename),status='unknown')
c     .   '/'//trim(filename),status='unknown', access='append')
         xq = xq_initial
         do i=1,it_max
          write(21,*)xq,PKPlus(i),err_Plus(i)
          xq = xq + xincr
         enddo
         close(21)
123        continue
       end
c ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
c      subroutine intPK_Delta(PKDel,err_Del)
c      implicit double precision (a-h,o-z)
c      dimension c(1:2)
c      character*50 name,mode,mode1
c      character*100 command,run_tag,dir_path,filename,filename1
c      common/energy/s
c      external flo2_Plus,flo2_PKDel,flo2_PKReg
c      common/usedalpha/AL,ge   
c      common/distribution/xq
c      dimension PKPlus(1:50),err_Plus(1:50)
c      dimension PKReg(1:50),err_Reg(1:50)
c      dimension PKDel(1:50),err_Del(1:50)
c      dimension PK(1:50),err(1:50)
c      
c
c      !input data card
c      open(unit=10,file='../run.vegas.dat',status='unknown')
c      do i=1,6
c      read (10,*)
c      enddo
c      read (10,*) pt1           ! vegas points     
c      read (10,*) its1          ! vegas iterations 
c      npt1 = pt1/10d0
c      close(10)
c      
c
c      open(unit=10,file='../param_card.dat',status='unknown')    
c      read (10,*) ge          ! [ 1/Alpha_ew ]
c      close(10)
c
c      open(unit=15,file='../run.machine.dat',status='unknown')
c      read (15,*) mid           ! machine id Tevatron:0 LHC:1
c      read (15,*) ecm           ! ecm
c      read (15,*) name          ! lhapdf set
c      read (15,*) it_max        ! it_max no of q for distribution
c      read (15,*) xq            ! initialise xq value
c      read (15,*) xincr         ! increment in Gev from xq 
c      read (15,*) run_tag
c      read (15,*) iprint        ! to print data in file
c      close(15)
c
cc ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        
c
c      open(unit=20,file='../output_files.dat',status='unknown')
c      do i =1,7
c      read (20,*) 
c      enddo
c      read (20,*) filename
c      close(20)
c      if (iprint .eq. 1) call output(run_tag,filename)            
cc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        
c
c
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[P and K terms from Here ] 
c        xq_initial = xq
c      call initpdfsetbyname(name)
c      Call initPDF(0)
c       s=ecm*ecm
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ delta functions ]
c        xq = xq_initial
c
c      mode = "Delta Functions"
c      call printframe0(mode)
c      call printframe1(pt1,its1)
c
c        do j=1,it_max
c
c      call printframe2(xq)
c
cc     -------------------------------------------------
c        call brm48i(40,0,0) 
c        call vsup(2,npt1,its1,flo2_PKDel,ai_lo2,sd,chi2)
cc     -------------------------------------------------
c
c          PKDel(j) = ai_lo2
c          err_Del(j)=sd
c
c      mode = "Delta"
c      call printframe3(mode,ai_lo2,sd,chi2)
c
c           xq=xq + xincr
c        enddo
c        xq = xq_initial
c
c       call printframe4(mode)
c        do j=1,it_max
c          write(*,'(i7,3e27.15,3e27.15)')
c     .             int(xq),PKDel(j),err_Del(j)
c          xq = xq + xincr
c        enddo
c
cc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[  * END * ]      
c
c
ccc ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
cc        if (iprint .ne. 1 ) goto 123
cc       open(unit=21,file='../summary/'//trim(run_tag)//
cc     .   '/'//trim(filename),status='unknown')
ccc     .   '/'//trim(filename),status='unknown', access='append')
cc         xq = xq_initial
cc         do i=1,it_max
cc          write(21,*)xq,PKDel(i),err_Del(i)
cc          xq = xq + xincr
cc         enddo
cc         close(21)
c123        continue
c       end
cc ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
