c~~~~~~~[	NIRAJ FILES OPEN HERE	   ]~~~~~~~~~~~~~~~O
	open(unit=560,file='nfile.dat',status='unknown')
	read(560,*)run_tag
	read(560,*)iprint
	close(560)
       if (iprint .eq. 1) then
        if (norder .eq. 0 ) then
         call output(run_tag,fname)
        open(unit=561,file='summary/'//trim(run_tag)//'/'
     .   //trim(fname),status='unknown')
        elseif (norder.eq.1) then
       call output(run_tag,fname1)
       open(unit=571,file='summary/'//trim(run_tag)//'/'
     .  //trim(fname1),status='unknown')
       call output(run_tag,fname2)
       open(unit=581,file='summary/'//trim(run_tag)//'/'
     .  //trim(fname2),status='unknown')
       call output(run_tag,fname3)
       open(unit=591,file='summary/'//trim(run_tag)//'/'
     .  //trim(fname3),status='unknown')
       call output(run_tag,fname4)
       open(unit=801,file='summary/'//trim(run_tag)//'/'
     .  //trim(fname4),status='unknown')
       call output(run_tag,fname5)
       open(unit=851,file='summary/'//trim(run_tag)//'/'
     .  //trim(fname5),status='unknown')
       endif
       endif


