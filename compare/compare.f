        program compare
        implicit double precision (a-h,o-z)
        dimension xqLO(1:50),xintLO(1:50),xqvir(1:50),xintvir(1:50),
     .    xqreal(1:50),xintreal(1:50),xqPKterm1(1:50),xintPKterm1(1:50)
     .                                ,xqPKterm2(1:50),xintPKterm2(1:50)

        
c       Leading Order

        open(unit=17,file='LO.dat',status='unknown')
        do i=1,20
        read(17,*) xqLO(i),xintLO(i)
        enddo
        close(17)
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
c     . achar(27)//'[0m'
        do i=1,20
c        write(*,'(i7,3e27.15)')int(xqLO(i)),xintLO(i)
        enddo

c       Virtual contribution

        open(unit=17,file='vir.dat',status='unknown')
        do i=1,20
        read(17,*) xqVir(i),xintVir(i)
        enddo
        close(17)
c        print*," "
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
c     . achar(27)//'[0m'
        do i=1,20
c          write(*,'(i7,3e27.15)')int(xqVir(i)),xintVir(i)
        enddo

c Ratio will be
c        print*," "
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","          ratio 
c     .virtual/LO",achar(27)//'[0m'
        do i=1,20
c          write(*,'(i7,3e27.15)')int(xqVir(i)),xintVir(i)/xintLO(i)
        enddo

c       Real emission Contribution        

        print*,"Real - Dipole"
        open(unit=17,file='real.dat',status='unknown')
        do i=1,20
        read(17,*) xqreal(i),xintreal(i)
        enddo
        close(17)

        print*," "
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
c     . achar(27)//'[0m'
        do i=1,20
c          write(*,'(i7,3e27.15)')int(xqreal(i)),xintreal(i)
        enddo

c        print*,"PK term1"
        open(unit=17,file='PKterm1.dat',status='unknown')
        do i=1,20
        read(17,*) xqPKterm1(i),xintPKterm1(i)
        enddo
        close(17)
c        print*," "
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
c     . achar(27)//'[0m'
        do i=1,20
c          write(*,'(i7,3e27.15)')int(xqPKterm1(i)),xintPKterm1(i)
        enddo

c        print*,"PK term2"
c        print*," "
        open(unit=17,file='PKterm2.dat',status='unknown')
        do i=1,20
        read(17,*) xqPKterm2(i),xintPKterm2(i)
        enddo
        close(17)
c        print*," "
c       write(*,*)achar(27)//'[1;32m'//"   xq"," ","           Integral",
c     . achar(27)//'[0m'
        do i=1,20
c          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),xintPKterm2(i)
        enddo


c       Total sig_NLO will be

        print*," "
       write(*,*)achar(27)//'[1;32m'//"   xq"," ","          Sigma_NLO",
     . achar(27)//'[0m'
        do i=1,20
          write(*,'(i7,3e27.15)')int(xqPKterm2(i)),xintPKterm2(i)+
     .  xintPKterm1(i)+xintvir(i)+xintreal(i)
        enddo
        end
