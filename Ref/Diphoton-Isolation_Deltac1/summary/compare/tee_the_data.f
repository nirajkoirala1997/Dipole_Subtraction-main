      program handle_output
      implicit double precision (a-h,o-z)
      character* 25 run_tag
      character* 50 mode

       call get_command_argument(1,mode)


      open(unit=15,file='../../nfile.dat',status='unknown')
      read(15,*)run_tag
      close(15)
      
      call system("mv ../output.dat ../"//trim(run_tag))
      end


