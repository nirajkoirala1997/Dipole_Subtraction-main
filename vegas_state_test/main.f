
      program eE2uUg_dipoleSubtraction 
      implicit double precision (a-h,o-z)
      common/recover/ireset

      INTEGER  seed
      LOGICAL state_exists,integral_state
      external fnlo3

      !input data card
      open(unit=10,file='run.vegas.dat',status='unknown')    
      read (10,*) pt1      
      npt1 = pt1
      read (10,*) its1    
      close(10)

! Check if state file exists
      state_exists = .FALSE.
c      INQUIRE(FILE='state.dat', EXIST=state_exists)
c      INQUIRE(FILE='integral.state', EXIST=integral_state)
      
      IF (state_exists) THEN
          ! Restore state from file
          OPEN(UNIT=10, FILE='state.dat', STATUS='OLD')
          READ(10, *) IJKLUT, NTOTUT, NTOT2T
          CLOSE(10)
          OPEN(UNIT=11, FILE='integral.state', STATUS='OLD')
          READ(11, *) IT,TI,TSI,AVGI,SD,CHI2A 
          CLOSE(11)
          CALL BRM48I(IJKLUT, NTOTUT, NTOT2T)
          ireset =1
      ELSE
          ! Initialize RNG with a new seed
          seed = 40 
          CALL BRM48I(seed, 0, 0)
c          call brm48i(40,0,0) ! initialize random number generator
          call vsup(1,npt1,its1,fnlo3,answer,sd,chi2)
      END IF
        write(*,*)'The answer is =', answer
        write(*,*)"Integral      =",answer,"+-",sd
        write(*,*)"with chisq    =",chi2

c----------------------------------------------------------------
c       call system('rm integral.state state.dat')
       end

      function fnlo3(xx,vwgt)
      implicit double precision (a-h,o-z)
      dimension xx(1)

	fnlo3 = 0d0
	delta = 1d-15
	if (xx(1) .le. 1d0-delta) then
      fnlo3 = 1d0+xx(1)**2
      fnlo3 = fnlo3/(1d0-xx(1))-2d0*dLog(1d0-xx(1))
	endif
     
      return
      end
