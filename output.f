        program output
        implicit double precision (a-h,o-z)
      open(unit=15,file='summary/summary2.dat',status='unknown')
      read (15,*) real_dip        !iorder no of q for distribution
      read (15,*) virt            ! initialise xq value
      read (15,*) PK            ! initialise xq value
      close(15)
      print*,"total=",real_dip+virt+PK
      end

