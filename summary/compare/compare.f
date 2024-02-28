        program compare
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPKterm1(1:50),xintPKterm1(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)
        character*100 run_tag
        character*50 name

        
c       Leading Order
      open(unit=15,file='../../run.machine.dat',status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable 
      read (15,*) run_tag               ! name of run directory to save output
      close(15)

      open(unit=15,file='../'//trim(run_tag)//'/run.machine.dat'
     . ,status='unknown')
      read (15,*) mid                   ! machine id Tevatron:0 LHC:1
      read (15,*) ecm                   ! ecm
      read (15,*) name                  !lhapdf set
      read (15,*) it_max                !it_max no of q for distribution
      read (15,*) xq_initial            ! initialise xq value
      read (15,*) step_size             ! size in the multiplle of loop variable
      read (15,*) run_tag               ! name of run directory to save output
      close(15)

      print*,"Reading Data from "//trim(run_tag)
      call sleep(1)


        open(unit=17,file='../'//trim(run_tag)//'/LO.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqLO(i),xintLO(i)
        enddo
        close(17)
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","        Integral_LO",
     . achar(27)//'[0m'
        do i=1,it_max
        write(*,'(i7,3e27.15)')int(xqLO(i)),xintLO(i)
        enddo

c       Virtual contribution

        open(unit=17,file='../'//trim(run_tag)//'/virtual.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqVir(i),xintVir(i)
c        read(17,*) xqVir(i),xintVir(i)
        enddo
        close(17)
        print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","       Integral_VIR",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqVir(i)),xintVir(i)
        enddo

c Ratio will be
c        print*," "
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","          ratio 
c     . LO_CH/LO","LO/LO_CH",achar(27)//'[0m'
c        do i=1,it_max
c       write(*,'(i7,3e27.15,3e27.15)')int(xqVir(i)),xintLO(i)/xintVir(i)
c     .             ,xintVir(i)/xintLO(i)
c        enddo
c       Real emission Contribution        

c        print*,"Real - Dipole"
        open(unit=17,file='../'//trim(run_tag)//'/real.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqreal(i),xintreal(i)
c        read(17,*) xqreal(i),xintreal(i)
        enddo
        close(17)

        print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","      Integral Real",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqreal(i)),xintreal(i)
        enddo

c ~~~~~~~~~~~~~~~~~~~~~~~~`PK terms

        open(unit=17,file='../'//trim(run_tag)//'/PK1.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqPKterm1(i),xintPKterm1(i)
c        read(17,*) xqPKterm1(i),xintPKterm1(i)
        enddo
        close(17)
        print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","      Integral PK1",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqPKterm1(i)),xintPKterm1(i)
        enddo

        open(unit=17,file='../'//trim(run_tag)//'/PK2.dat',
     .     status='unknown')
        do i=1,it_max
        read(17,*) xqPKterm2(i),xintPKterm2(i)
c        read(17,*) xqPKterm2(i),xintPKterm2(i)
        enddo
        close(17)
        print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","       Integral PK2",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),xintPKterm2(i)
        enddo

c ~~~~~~~~~~~~~~~~~~~~~~~~`PK terms end


c        print*," "
c        print*,"Sigma_Chinmoy"
c        print*," "
c        open(unit=17,file='LO_2chinmoy.dat',status='unknown')
c        do i=1,it_max
c        read(17,*) xqLO(i),xintLO(i)
c        enddo
c        close(17)
c        print*," "
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","       ChinIntegral",
c     . achar(27)//'[0m'
c        do i=1,it_max
c          write(*,'(i7,3e27.15)')int(xqLO(i)),xintLO(i)
c        enddo


c       Total sig_NLO will be

120     print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","          Sigma_NLO",
     . achar(27)//'[0m'
        do i=1,it_max
          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),xintPKterm2(i)+
     .  xintPKterm1(i)+xintvir(i)+xintreal(i)!+xintLO(i)
c          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),
c     .  xintvir(i)+xintreal(i)+xintLO(i)
        enddo
        end
