      program handle_output
       implicit double precision (a-h,o-z)
       character*50 name,filename,run_tag,filename_tmp,mode,tag
     .  ,filename1,filename2,filename3,filename4

       open(unit=15,file='../../run.machine.dat',status='unknown')
       read (15,*) mid           ! machine id Tevatron:0 LHC:1
       read (15,*) ecm           ! ecm
       read (15,*) name          ! lhapdf set
       read (15,*) it_max        ! it_max no of q for distribution
       read (15,*) xq            ! initialise xq value
       read (15,*) xincr         ! increment in Gev from xq 
       read (15,*) run_tag
       read (15,*) iprint        ! to print data in file
       close(15)

       call get_command_argument(1,mode)
       call get_command_argument(2,tag)
 

       open(unit=20,file='../../output_files.dat',status='unknown')
       read (20,*) filename1
       read (20,*) filename2
       read (20,*) filename3
       read (20,*) filename4
       close(20)

       if (tag .eq. 'dipole') filename = filename1
       if (tag .eq. 'virtual') filename = filename2
       if (tag .eq. 'PK') filename = filename3
       if (tag .eq. 'LO') filename = filename4
       if (iprint .eq. 1 ) then
       call system("cd ../../trash/broken && cat "//trim(mode)// " >>
     .  ../../summary/"//trim(run_tag)//"/"//trim(filename)//"
     .   && rm "//trim(mode))
        filename_tmp = "../"//trim(run_tag)//"/"//trim(filename)

       call system(
     . 'sed -i -r "s/\\x1B\\[[0-9;]*[a-zA-Z]//g" '//trim(filename_tmp))
       endif
       
       end

