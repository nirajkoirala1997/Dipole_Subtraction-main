c       PROGRAM MAIN
c       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c       INTEGER ND
c       XMH=300D0 
c       XMS=4000D0 
c       ND=3
c       CALL FACTGR(XMH,XMS,ND,FGR)
c       WRITE(*,*)"FGR=",FGR
c       END PROGRAM MAIN
C  FGR=   8.9224175345902017E-015
C FLUX Prefactors are defined here
C******************************************************

       SUBROUTINE XmFACTUD(XMH,FUUN,FDDN)
*      --------------------------------------------------
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON/VMASS/XMZ,XMW,XMMH,GMZ,GMW
       COMMON/WENBERG/SW2,CW2,ALEM
       COMMON/AHIGGS/NHIGGS
       COMMON/PROCESS/NPROC
       DOUBLE PRECISION HIGGS
       EXTERNAL HIGGS
C
       XMH2=XMH*XMH
       XMZ2=XMZ*XMZ
       GMZ2=GMZ*GMZ
       ALEM2=ALEM*ALEM
       XMMH2=XMMH*XMMH
C
       QE=-1.0D0
       QU=2.0D0/3.0D0
       QD=-1.0D0/3.0D0
C
       GEVZ=-1.0D0/4.0D0-SW2*QE
       GUVZ=1.0D0/4.0D0-SW2*QU
       GDVZ=-1.0D0/4.0D0-SW2*QD
C
       GEAZ=1.0D0/4.0D0
       GUAZ=-1.0D0/4.0D0
       GDAZ=1.0D0/4.0D0
C
       TU=1.0d0/2.0d0
       TD=-1.0d0/2.0d0
       
       IF (NHIGGS.EQ.1) THEN

       FACTORZ = (3.0d0*XMZ2/XMH**2)*(HIGGS(XMZ2,XMMH2,XMH**2))**0.5d0*
     &       (1.0d0+1.0d0/12.0d0*HIGGS(XMZ2,XMMH2,XMH**2)*XMH**2/XMZ2)/
     &       (GEAZ*GEAZ+GEVZ*GEVZ)/4.0d0

       ELSEIF (NHIGGS.EQ.0) THEN

       FACTORZ = 1.0D0

       ENDIF

       IF (NPROC.EQ.0) THEN

       T1UN=QU*(2.0D0*QU+3.0D0*QD)
       T21UN=-XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         QU*GEVZ*(-1.0D0/4.0D0-SW2*(2.0D0*QU+3.0D0*QD))
       T22UN=-XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         GEVZ*GUVZ*(2.0D0*QU+3.0D0*QD)
       T3UN=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(-1.0D0/8.0D0+1.0D0/12.0D0*SW2
     &                           +2.0D0/9.0D0*SW2**2)
C
       T1DN=QD*(2.0D0*QU+3.0D0*QD)
       T21DN=-XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         QD*GEVZ*(-1.0D0/4.0D0-SW2*(2.0D0*QU+3.0D0*QD))
       T22DN=-XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         GEVZ*GDVZ*(2.0D0*QU+3.0D0*QD)
       T3DN=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(1.0D0/8.0D0-1.0D0/9.0D0*SW2**2)
C
*test
       FUUN=0.0d0
       FDDN=0.0d0
*test
       FUUN=FUUN+4.0D0/3.0D0*ALEM2/XMH2*(T1UN+T21UN+T22UN+T3UN)
       FDDN=FDDN+4.0D0/3.0D0*ALEM2/XMH2*(T1DN+T21DN+T22DN+T3DN)

       ELSEIF(NPROC.EQ.1) THEN

       ZBOSON = FACTORZ*((4.0D0/3.0D0*ALEM2/XMH**2)*
     &       (XMH**4.0d0/((XMH**2 - XMZ2)**2+XMZ2*GMZ2)))

        TU3N = 1.0d0/CW2/CW2/SW2/SW2*(GEVZ*GEVZ+GEAZ*GEAZ)*
     &        (-1.0D0/8.0D0+1.0D0/12.0D0*SW2+2.0D0/9.0D0*SW2**2)

        T3DN = 1.0d0/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(1.0D0/8.0D0-1.0D0/9.0D0*SW2**2)

       FUUN=(T3UN)*ZBOSON
       FDDN=(T3DN)*ZBOSON
        
c       write(*,*)'FUUN, FDDN = ',FUUN,FDDN

       ENDIF


       RETURN
       END

       SUBROUTINE FACTUD(XMH,FUU,FDD,FUDV,FUDA,FDUV,FDUA)
*      --------------------------------------------------
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON/VMASS/XMZ,XMW,XMMH,GMZ,GMW
       COMMON/WENBERG/SW2,CW2,ALEM
       COMMON/PROCESS/NPROC
       COMMON/AHIGGS/NHIGGS
       DOUBLE PRECISION HIGGS
       EXTERNAL HIGGS

c       write(*,*)'XMZ, XMW, GMZ, ALEM, XMMH',XMZ,XMW,GMZ,ALEM,XMMH
C
       XMH2=XMH*XMH
       XMZ2=XMZ*XMZ
       XMW2=XMW*XMW
       GMZ2=GMZ*GMZ
       GMW2=GMW*GMW
       ALEM2=ALEM*ALEM
       XMMH2=XMMH*XMMH
C
       QE=-1.0D0
       QU=2.0D0/3.0D0
       QD=-1.0D0/3.0D0
       SR2 = DSQRT(2.0D0)
C
       GEVZ=-1.0D0/4.0D0-SW2*QE
       GUVZ=1.0D0/4.0D0-SW2*QU
       GDVZ=-1.0D0/4.0D0-SW2*QD

       GEAZ=1.0D0/4.0D0
       GUAZ=-1.0D0/4.0D0
       GDAZ=1.0D0/4.0D0

       GEVW=1.0d0/2.0d0/SR2
       GUVW=1.0d0/2.0d0/SR2
       GDVW=1.0d0/2.0d0/SR2
C
       GEAW=1.0d0/2.0d0/SR2
       GUAW=1.0d0/2.0d0/SR2
       GDAW=1.0d0/2.0d0/SR2

       IF (NHIGGS.EQ.1) THEN

c       FACTOR = (3.0d0*XMW2/XMH**2)*(HIGGS(XMW2,XMMH2,XMH**2))**0.5d0*
c     &(1.0d0+1.0d0/12.0d0*HIGGS(XMW2,XMMH2,XMH**2)*XMH**2/XMW2)/CW2


       FACTORZ = (3.0d0*XMZ2/XMH**2)*(HIGGS(XMZ2,XMMH2,XMH**2))**0.5d0*
     &       (1.0d0+1.0d0/12.0d0*HIGGS(XMZ2,XMMH2,XMH**2)*XMH**2/XMZ2)/
     &       (GEAZ*GEAZ+GEVZ*GEVZ)/4.0d0

       FACTORW = (3.0d0*XMW2/XMH**2)*(HIGGS(XMW2,XMMH2,XMH**2))**0.5d0*
     &       (1.0d0+1.0d0/12.0d0*HIGGS(XMW2,XMMH2,XMH**2)*XMH**2/XMW2)/
     &       (GEAW*GEAW+GEVW*GEVW)/4.0d0

       ELSEIF (NHIGGS .EQ. 0) THEN

       FACTORZ = 1.0d0
       FACTORW = 1.0d0

       ENDIF

       IF (NPROC.EQ.0) THEN

       T1U=QU*QU
       T2U=-2.0D0*XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         QU*GEVZ*GUVZ
       T3U=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GUVZ*GUVZ+GUAZ*GUAZ)
C
       T1UD=QU*QD
       T2UD=-2.0D0*XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         GEVZ*(QU*GDVZ+QD*GUVZ)/2.0D0
       T3UDV=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GUVZ*GDVZ)
       T3UDA=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GUAZ*GDAZ)
C
       T1D=QD*QD
       T2D=-2.0D0*XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         QD*GEVZ*GDVZ
       T3D=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GDVZ*GDVZ+GDAZ*GDAZ)
C
       T1DU=QD*QU
       T2DU=-2.0D0*XMH2*(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &         GEVZ*(QD*GUVZ+QU*GDVZ)/2.0D0
       T3DUV=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GDVZ*GUVZ)
       T3DUA=XMH2*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GDAZ*GUAZ)
C
       S1U=(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &      QU*QE*GUAZ*GEAZ
       S2U=2.0D0*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2/CW2/SW2*
     &      GUVZ*GEVZ*GUAZ*GEAZ
       S1D=(XMH2-XMZ2)/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2*
     &      QD*QE*GDAZ*GEAZ
       S2D=2.0D0*XMH2/((XMH2-XMZ2)**2+XMZ2*GMZ2)/CW2/SW2/CW2/SW2*
     &      GDVZ*GEVZ*GDAZ*GEAZ

C       
c       FSMU=4.0D0/3.0D0*ALEM2/XMH2*(T1U+T2U+T3U)
c       FSMD=4.0D0/3.0D0*ALEM2/XMH2*(T1D+T2D+T3D)
C
       FUU=4.0D0/3.0D0*ALEM2/XMH2*(T1U+T2U+T3U)
       FDD=4.0D0/3.0D0*ALEM2/XMH2*(T1D+T2D+T3D)
C
       FUDV=4.0D0/3.0D0*ALEM2/XMH2*(T1UD+T2UD+T3UDV)
       FUDA=4.0D0/3.0D0*ALEM2/XMH2*(T1UD+T2UD+T3UDA)
C
       FDUV=4.0D0/3.0D0*ALEM2/XMH2*(T1DU+T2DU+T3DUV)
       FDUA=4.0D0/3.0D0*ALEM2/XMH2*(T1DU+T2DU+T3DUA)

       ELSEIF(NPROC.EQ.1)THEN

       ZBOSON = FACTORZ*((4.0D0/3.0D0*ALEM2/XMH**2)*
     &       (XMH**4.0d0/((XMH**2 - XMZ2)**2+XMZ2*GMZ2)))

       T3U = 1.0D0/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GUVZ*GUVZ+GUAZ*GUAZ)
C
       T3D = 1.0D0/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GDVZ*GDVZ+GDAZ*GDAZ)

       T3UDV=1.0d0/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GUVZ*GDVZ)

       T3UDA=1.0d0/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GUAZ*GDAZ)

       T3DUV=1.0d0/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GDVZ*GUVZ)

       T3DUA=1.0d0/CW2/CW2/SW2/SW2*
     &        (GEVZ*GEVZ+GEAZ*GEAZ)*(GDAZ*GUAZ)

C
       FUU = (T3U)*ZBOSON
       FDD = (T3D)*ZBOSON
c       write(*,*)'FUU, FDD =', T3U, T3D

       FUDV=(T3UDV)*ZBOSON
       FUDA=(T3UDA)*ZBOSON
C
       FDUV=(T3DUV)*ZBOSON
       FDUA=(T3DUA)*ZBOSON

       ELSEIF(NPROC.EQ.2.OR.NPROC.EQ.3.OR.NPROC.EQ.4)THEN

       WBOSON = FACTORW*((4.0D0/3.0D0*ALEM2/XMH**2)*
     &       (XMH**4.0d0/((XMH**2 - XMW2)**2+XMW2*GMW2)))

       T3U = 1.0D0/SW2/SW2*
     &        (GEVW*GEVW+GEAW*GEAW)*(GUVW*GUVW+GUAW*GUAW)
C
       T3D = 1.0D0/SW2/SW2*
     &        (GEVW*GEVW+GEAW*GEAW)*(GDVW*GDVW+GDAW*GDAW)

       T3UDV=1.0d0/SW2/SW2*
     &        (GEVW*GEVW+GEAW*GEAW)*(GUVW*GDVW)

       T3UDA=1.0d0/SW2/SW2*
     &        (GEVW*GEVW+GEAW*GEAW)*(GUAW*GDAW)

       T3DUV=1.0d0/SW2/SW2*
     &        (GEVW*GEVW+GEAW*GEAW)*(GDVW*GUVW)

       T3DUA=1.0d0/SW2/SW2*
     &        (GEVW*GEVW+GEAW*GEAW)*(GDAW*GUAW)

C
       FUU = (T3U)*WBOSON
       FDD = (T3D)*WBOSON
c       write(*,*)'FUU, FDD =', T3U, T3D

       FUDV=(T3UDV)*WBOSON
       FUDA=(T3UDA)*WBOSON
C
       FDUV=(T3DUV)*WBOSON
       FDUA=(T3DUA)*WBOSON

       ENDIF
C
       RETURN
       END

       REAL*8 FUNCTION HIGGS(X,Y,Z)
*      ------------------------
       IMPLICIT REAL*8 (A-H,O-Z)
C
       HIGGS=( 1 - X/Z - Y/Z )**2.0d0 - 4.0d0*X*Y/Z/Z
C
       RETURN
       END



**********************************************************************
************************|--------------------|************************
************************| Program: DiGamma.f |************************
************************|--------------------|************************
**********************************************************************
*****************|-------------------------------|********************
*****************| Last modified:  June 1, 2000  |********************
*****************|-------------------------------|********************
**********************************************************************
********|----------------------------------------------------|********
********| This function calculates the Euler Digamma function|********
********| for a real argument using the series expansion:    |********
********|                                                    |********
********|  \psi(x) = - \gamma                                |********
********|             + \sum_{n=1}^{\infty} (x-1)/( n(n+x-1) |********
********|                                                    |********
********| given by Abramowitz and Stegun.                    |********
********| \gamma is Euler's constant = 0.57721566            |********
********|----------------------------------------------------|********
********| Warning: This does not converge for x = 0,-1,-2,.. |********
********|----------------------------------------------------|********
**********************************************************************
*--------------------------------------------------------------------*
          Function DiGamma (X)
*--------------------------------------------------------------------*
*  Declarations:
*  ~~~~~~~~~~~~
          Implicit None
          Integer n
          Real*8  X, DiGamma, EulConst, Term, Eps
          Data EulConst/0.57721566d0/
*--------------------------------------------------------------------*
             n = 0
                             DiGamma = - EulConst
 99          n = n +1
          Term = (X - 1.d0) /n/(n - 1.d0 + X)
          Eps  = DAbs(Term/DiGamma)
          If( Eps. Gt. 1.D-4) Then
                             DiGamma = DiGamma + Term
                             Go To 99
          EndIf
*--------------------------------------------------------------------*
          Return
          End
**********************************************************************
