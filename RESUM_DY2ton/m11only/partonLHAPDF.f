CC	 Evolves LHAPDF pdf sets for different Q and x through evolvepdf function of LHAPDF (23.09.2017)
CC       Note that down has PDF-1, up PDF-2

      SUBROUTINE STRUCT(X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,
     &       STRB,CHM,CHMB,BOT,BOTB,GLU)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 df(-6:6),ff(-6:6)

      call evolvepdf(X,SCALE,ff)

      glu  = ff(0)
      dn   = ff(1)
      dsea = ff(-1)
      up   = ff(2)
      usea = ff(-2)
      str  = ff(3)
      chm  = ff(4)
      bot  = ff(5)
      upv  = up-usea
      dnv  = dn-dsea
      strb = ff(-3)
      chmb  = ff(-4)
      botb  = ff(-5)

      RETURN
      END

