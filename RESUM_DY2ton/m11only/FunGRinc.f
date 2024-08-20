*----------------------------------------------------------------------------------------
       DOUBLE PRECISION FUNCTION REG(YY) !! this is the functions after the N integration
*      ----------------------------
       IMPLICIT NONE

       DATA PI/3.141592653589793238462643D0/

       DOUBLE PRECISION YY
       DIMENSION YY(5)


      DOUBLE PRECISION TAUH,xINT1,xINT2,PI
      DOUBLE PRECISION Y1,Y2,Y3,Y4,Y5,X1,X2,X3,X4,UL1,PHI1
      DOUBLE COMPLEX ONE,ZERO,II,XIM,C1,EXPIPHI,xN1
      DOUBLE PRECISION XH,X1P,X2P,X1M,X2M
      DOUBLE PRECISION UV1,DV1,US1,DS1,ST1,CH1,BOT1,GL1,ST1B,CH1B,BOT1B
      DOUBLE PRECISION UV2,DV2,US2,DS2,ST2,CH2,BOT2,GL2,ST2B,CH2B,BOT2B
      DOUBLE PRECISION UV1P,DV1P,US1P,DS1P,ST1P,CH1P,BOT1P,GL1P,ST1BP
     &                 ,CH1BP,BOT1BP 
      DOUBLE PRECISION UV2P,DV2P,US2P,DS2P,ST2P,CH2P,BOT2P,GL2P,ST2BP
     &                 ,CH2BP,BOT2BP 
      DOUBLE PRECISION UV1M,DV1M,US1M,DS1M,ST1M,CH1M,BOT1M,GL1M,ST1BM
     &                 ,CH1BM,BOT1BM 
      DOUBLE PRECISION UV2M,DV2M,US2M,DS2M,ST2M,CH2M,BOT2M,GL2M,ST2BM 
     &                 ,CH2BM,BOT2BM 
      DOUBLE PRECISION U1,U1B,D1,D1B,U2,U2B,D2,D2B 
      DOUBLE PRECISION U1P,U1BP,D1P,D1BP,U2P,U2BP,D2P,D2BP 
      DOUBLE PRECISION U1M,U1BM,D1M,D1BM,U2M,U2BM,D2M,D2BM 
      DOUBLE PRECISION FU1,FU1B,FD1,FD1B,FU2,FU2B,FD2,FD2B
      DOUBLE PRECISION FCH1,FCH2,FST1,FST2,FBOT1,FBOT2
     &                 ,FST1B,FST2B,FCH1B,FCH2B,FBOT1B,FBOT2B
      DOUBLE PRECISION CCH1,CCH2,SST1,SST2,BBOT1,BBOT2,SST1B,SST2B
     &                 ,CCH1B,CCH2B,BBOT1B,BBOT2B
      DOUBLE PRECISION UU1,UU1B,DD1,DD1B,UU2,UU2B,DD2,DD2B 
      DOUBLE PRECISION FGL1,FGL2,GGL1,GGL2

      DOUBLE PRECISION FLUX,HQQBd,HQQBu 
      DOUBLE PRECISION BORN,FLO2,VIRT,VIRT1,VIRT2,VIRTN
      DOUBLE PRECISION Qmax,Qmin,Q2,XJAC,XJAC1,XJAC2,XJAC3,XEPS 
      DOUBLE COMPLEX FUNC

      DOUBLE PRECISION VUD, VUS, VUB, VCD, VCS, VCB
      DOUBLE PRECISION VUD2, VUS2, VUB2, VCD2, VCS2, VCB2
      DOUBLE PRECISION HELP1, HELP2

      DOUBLE PRECISION GF,XMZ,XMW,SW2,CW2,GMZ,GMW,ALEM,XMMH,stw,ctw    
      real*8 gau,gvu,gad,gvd
      DOUBLE PRECISION g01zz  

!!! Common block
      INTEGER MODL,MODE,MODP,CHANNEL,NF,ORD,INFV,NPROC,NHIGGS
      DOUBLE PRECISION RS,CA,CF,TF,XNC,XMT,XMB,XMC,Q,XMUR,XMUF,AS
      DOUBLE PRECISION XLFR,XLQF,XLQR
      DOUBLE PRECISION XLAMBDA,XKAP
      CHARACTER*50 PREFIX

      COMMON/ADDRS/MODL
      COMMON/ISET/MODE,MODP
      COMMON/MACHINE/RS
      COMMON/CHANNEL/CHANNEL
      COMMON/COLORS/CA,CF,TF,XNC
      COMMON/FLAVOUR/NF
      COMMON/PDFGROUP/PREFIX
      COMMON/XMASS/XMT,XMB,XMC
      COMMON/NORDER/ORD
      COMMON/FLAVOUR/INFV
      COMMON/XXMH/Q,XMUR,XMUF
      COMMON/ALS/AS
      COMMON/COUPLING/XLAMBDA,XKAP
      COMMON/SCALES/XLFR,XLQF,XLQR
      COMMON/CKM/VUD,VUS,VUB,VCD,VCS,VCB
      COMMON/PROCESS/NPROC
      COMMON/AHIGGS/NHIGGS
      COMMON/Qvalue/Qmin,Qmax,XEPS
      COMMON/WENBERG/SW2,CW2,ALEM
      COMMON/gzzs/born,virt1,virt2
c      EXTERNAL FUNC
c      EXTERNAL FUNC,FLO2,VIRT
!!! ADD
       REAL*8 c00,ams,amh,aam1,c0,aamh
c      COMMON/CHANNEL/CHANNEL
c       COMMON/ADDRS/MODL
c       COMMON/ISET/MODE,MODP
c       common/xpar/c00,ams,amh
       COMMON/rs_par/aam1,c0,aamh
       DOUBLE PRECISION XMS,FGR
       INTEGER NOD
       COMMON/ADDPAR/XMS,NOD
       DOUBLE PRECISION FUU,FDD,FUDV,FUDA,FDUV,FDUA,FUUN,FDDN

      INTEGER IMODE
      COMMON/SVRES/IMODE

       Y1 = YY(1)
       Y2 = YY(2)
       Y3 = YY(3)
       Y4 = YY(4)
       Y5 = YY(5)


       XJAC3 = (Qmax-Qmin)
       Q = XJAC3*Y5+Qmin

       Q2 = Q*Q
       TAUH=Q2/RS/RS

       X1 = (1D0-TAUH)*Y1    + TAUH       
       XJAC1 = (1D0-TAUH)

       X2 = (1D0-TAUH/X1)*Y2 + TAUH/X1 
       XJAC2 = (1D0-TAUH/X1)

       X3 = TAUH/X1/X2               
       X4 = Y4

       XJAC = XJAC1*XJAC2*XJAC3

C      Born PS inegration
c       IF(ORD.EQ.0) THEN
       call Event(X1,TAUH/X1,X4,flo2,virt,virtN)          !! Numerical BORN in EXACT
       BORN = flo2                                !! Numerical BORN in EXACT
       BORN = BORN*2d0*Q/RS/RS        !! why Q?`

       virt1 = virt-1             !! virtual by born

       ONE  = DCMPLX(1.0D0,0.0D0)
       II   = (0.0D0,1.0D0)
       XIM  = DCMPLX(0.0D0,1.0D0)
       UL1  = 75.0D0
       C1   = DCMPLX(1.9D0,0.0D0)
c cos is zero
       PHI1 = 3.0D0*PI/4.0D0

       EXPIPHI = EXP(II*PHI1)
       xN1     = UL1*Y3*EXPIPHI+C1 
c       xINT1   = UL1*DIMAG(EXPIPHI/PI * EXP(-xN1*LOG(X3))* FUNC(xN1)) !??
       Call MellinInv(xN1,virt1,Func)
       xINT1   = UL1*DIMAG(EXPIPHI/PI * EXP(-xN1*LOG(X3))* FUNC) !??

       If (ORD.eq.2) then
       Call MellinInv(xN1,virt2,Func)
       xINT2   = UL1*DIMAG(EXPIPHI/PI * EXP(-xN1*LOG(X3))* FUNC) !??
       ENDIF

       XH  = 1.0D-7
       X1P = X1 + XH
       X1M = X1 - XH
       X2P = X2 + XH
       X2M = X2 - XH
       CALL STRUCT(X1,XMUF,MODE,UV1,DV1,US1,DS1,ST1,ST1B,CH1,
     &   CH1B,BOT1,BOT1B,GL1)
       CALL STRUCT(X2,XMUF,MODE,UV2,DV2,US2,DS2,ST2,ST2B,CH2,
     &   CH2B,BOT2,BOT2B,GL2)
       CALL STRUCT(X1P,XMUF,MODE,UV1P,DV1P,US1P,DS1P,ST1P,ST1BP,CH1P,
     &   CH1BP,BOT1P,BOT1BP,GL1P)
       CALL STRUCT(X2P,XMUF,MODE,UV2P,DV2P,US2P,DS2P,ST2P,ST2BP,CH2P,
     &   CH2BP,BOT2P,BOT2BP,GL2P)
       CALL STRUCT(X1M,XMUF,MODE,UV1M,DV1M,US1M,DS1M,ST1M,ST1BM,CH1M,
     &   CH1BM,BOT1M,BOT1BM,GL1M)
       CALL STRUCT(X2M,XMUF,MODE,UV2M,DV2M,US2M,DS2M,ST2M,ST2BM,CH2M,
     &   CH2BM,BOT2M,BOT2BM,GL2M)
*
C 1 =  PROTON, 2 = PROTON

       U1   = UV1+US1
       U1B  = US1
       D1   = DV1+DS1
       D1B  = DS1
       U2   = UV2+US2
       U2B  = US2
       D2   = DV2+DS2
       D2B  = DS2

       U1P  = UV1P+US1P
       U1BP = US1P
       D1P  = DV1P+DS1P
       D1BP = DS1P
       U2P  = UV2P+US2P
       U2BP = US2P
       D2P  = DV2P+DS2P
       D2BP = DS2P

       U1M  = UV1M+US1M
       U1BM = US1M
       D1M  = DV1M+DS1M
       D1BM = DS1M
       U2M  = UV2M+US2M
       U2BM = US2M
       D2M  = DV2M+DS2M
       D2BM = DS2M

c      useful formula
c      dy 
c      d/dx (x f(x)) = - int dN x^(-N)  (N-1) f(N)
c      d/dx ( x d/dx (x f(x)) ) =  int dN x^(-N)  (N-1)^2 f(N)
c      gluon
c      x d/dx (x f(x)) = - int dN x^(-N)  N f(N+1)
c      x d/dx ( x d/dx (x f(x)) ) =  int dN x^(-N)  N^2 f(N+1)


       FU1   = (U1P   - U1M)/2.0D0/XH
       FU1B  = (U1BP  - U1BM)/2.0D0/XH
       FU2   = (U2P   - U2M)/2.0D0/XH
       FU2B  = (U2BP  - U2BM)/2.0D0/XH
       FD1   = (D1P   - D1M)/2.0D0/XH
       FD1B  = (D1BP  - D1BM)/2.0D0/XH
       FD2   = (D2P   - D2M)/2.0D0/XH
       FD2B  = (D2BP  - D2BM)/2.0D0/XH
       FCH1  = (CH1P  - CH1M)/2.0D0/XH
       FCH1B  = (CH1BP  - CH1BM)/2.0D0/XH
       FCH2  = (CH2P  - CH2M)/2.0D0/XH
       FCH2B  = (CH2BP  - CH2BM)/2.0D0/XH
       FST1  = (ST1P  - ST1M)/2.0D0/XH
       FST1B  = (ST1BP  - ST1BM)/2.0D0/XH
       FST2  = (ST2P  - ST2M)/2.0D0/XH
       FST2B  = (ST2BP  - ST2BM)/2.0D0/XH
       FBOT1 = (BOT1P - BOT1M)/2.0D0/XH
       FBOT1B = (BOT1BP - BOT1BM)/2.0D0/XH
       FBOT2 = (BOT2P - BOT2M)/2.0D0/XH
       FBOT2B = (BOT2BP - BOT2BM)/2.0D0/XH

       FGL1  = (GL1P - GL1M)/2.0D0/XH
       FGL2  = (GL2P - GL2M)/2.0D0/XH

       UU1    = FU1   + X1*(U1P   + U1M   - 2.0D0*U1)/XH/XH
       UU1B   = FU1B  + X1*(U1BP  + U1BM  - 2.0D0*U1B)/XH/XH
       UU2    = FU2   + X2*(U2P   + U2M   - 2.0D0*U2)/XH/XH
       UU2B   = FU2B  + X2*(U2BP  + U2BM  - 2.0D0*U2B)/XH/XH
       DD1    = FD1   + X1*(D1P   + D1M   - 2.0D0*D1)/XH/XH
       DD1B   = FD1B  + X1*(D1BP  + D1BM  - 2.0D0*D1B)/XH/XH
       DD2    = FD2   + X2*(D2P   + D2M   - 2.0D0*D2)/XH/XH
       DD2B   = FD2B  + X2*(D2BP  + D2BM  - 2.0D0*D2B)/XH/XH
       CCH1   = FCH1  + X1*(CH1P  + CH1M  - 2.0D0*CH1)/XH/XH
       CCH2   = FCH2  + X2*(CH2P  + CH2M  - 2.0D0*CH2)/XH/XH
       SST1   = FST1  + X1*(ST1P  + ST1M  - 2.0D0*ST1)/XH/XH
       SST1B  = FST1B  + X1*(ST1BP  + ST1BM  - 2.0D0*ST1B)/XH/XH
       SST2   = FST2  + X2*(ST2P  + ST2M  - 2.0D0*ST2)/XH/XH
       SST2B  = FST2B  + X2*(ST2BP  + ST2BM  - 2.0D0*ST2B)/XH/XH
       BBOT1  = FBOT1 + X1*(BOT1P + BOT1M - 2.0D0*BOT1)/XH/XH
       BBOT2  = FBOT2 + X2*(BOT2P + BOT2M - 2.0D0*BOT2)/XH/XH

       GGL1  = FGL1  + X1*(GL1P  + GL1M  - 2.0D0*GL1)/XH/XH 
       GGL2  = FGL2  + X2*(GL2P  + GL2M  - 2.0D0*GL2)/XH/XH

        ctw = dsqrt(CW2)                            !Cos(thetaW)
        stw = dsqrt(sw2)

        gvu = 0.5d0*(0.5d0-2.0d0*sw2*2.0d0/3.0d0)/(stW*ctW)
        gvd = 0.5d0*(-0.5d0-2.0d0*sW2*(-1.0d0/3.0d0))/(stW*ctW)
        gau = 1.0d0/4d0/ctw/stw
        gad = gau

        ! FLUXES

c        FLUX = (  (FU1*FU2B+FU1B*FU2+2.0d0*FCH1*FCH2)*
c     &          (gau**4+6.0d0*gau**2*gvu**2+gvu**4)
c     &    +  (FD1*FD2B+FD1B*FD2+2.0d0*FST1*FST2+2.0d0*FBOT1*FBOT2)*
c     &     (gad**4+6.0d0*gad**2*gvd**2+gvd**4)   )

        FLUX = (  (FU1*FU2B+FU1B*FU2+FCH1*FCH2B+FCH1B*FCH2)*
     &          (gau**4+6.0d0*gau**2*gvu**2+gvu**4)
     &    + (FD1*FD2B+FD1B*FD2+FST1*FST2B+FST1B*FST2+FBOT1*FBOT2B
     &    + FBOT1B*FBOT2)*
     &     (gad**4+6.0d0*gad**2*gvd**2+gvd**4)   )

        REG = BORN*FLUX*xINT1*XJAC/2D0/XEPS/X1/X2

       RETURN

       END

*----------------------------------------------------------------------
**  The function to do Mellin inversion 
      Subroutine MellinInv(N,virt1,Func)
C     ---------------------------
      IMPLICIT NONE
!!! Common block
      INTEGER MODL,MODE,MODP,CHANNEL,NF,ORD,INFV
      DOUBLE PRECISION RS,CA,CF,TF,XNC,XMT,XMB,XMC,Q,XMUR,XMUF,AS
      DOUBLE PRECISION XLFR,XLQF,XLQR
      DOUBLE PRECISION XLAMBDA,XKAP
      CHARACTER*50 PREFIX

      COMMON/ADDRS/MODL
      COMMON/ISET/MODE,MODP
      COMMON/MACHINE/RS
      COMMON/CHANNEL/CHANNEL
      COMMON/COLORS/CA,CF,TF,XNC
      COMMON/FLAVOUR/NF
      COMMON/PDFGROUP/PREFIX
      COMMON/XMASS/XMT,XMB,XMC
      COMMON/NORDER/ORD
      COMMON/FLAVOUR/INFV
      COMMON/XXMH/Q,XMUR,XMUF
      COMMON/ALS/AS
      COMMON/COUPLING/XLAMBDA,XKAP
      COMMON/SCALES/XLFR,XLQF,XLQR
      INTEGER IMODE
      COMMON/SVRES/IMODE
c      COMMON/gzzs/born,virt1,virt2

      DOUBLE PRECISION bt0inv1,bt0inv2,bt0inv3,bt0inv4,bt0inv5
      DOUBLE PRECISION lfr,lqr,N4,nfv
      DOUBLE PRECISION CHg0tll,CHg0tnll,CHg0tnnll,CHg0tnnnll
      DOUBLE PRECISION CHg0tnnnll1,CHg0tnnnll2

      DOUBLE PRECISION BORN,FLO2,VIRT,VIRT1 
      DOUBLE PRECISION TAUH,x1,x4,g01zz

      DOUBLE COMPLEX N,NB,Nm1,NINV2,NINV4,ONE
      DOUBLE COMPLEX foll,fonll,fonnll,fonnnll
      DOUBLE COMPLEX fonnnll1,fonnnll2
      DOUBLE COMPLEX gll,gnll,gnnll,gnnnll
      DOUBLE COMPLEX lnNb,w,wm,winv,wminv,wminv2,zlwm,zlwm2
      DOUBLE COMPLEX ZINV,Func
      include 'constants.h'

      ONE   = DCMPLX(1.0D0,0.0D0)
      NB    = DEXP(GE)*N   !!! Nbar = Exp[GE]*N ! GE is EulerGamma
      Nm1   = N-ONE
         NINV2 = ONE/Nm1/Nm1
         NINV4 = ONE/Nm1/Nm1/Nm1/Nm1


CCCC   
      lnNb = ZLOG(NB)
      w       = 2D0*as*bt0*lnNb
      wm      = one - w
      winv    = 1.0D0/w
      wminv   = 1.0D0/wm
      wminv2  = 1.0D0/wm**2
      zlwm    = zlog(wm)
      zlwm2   = zlwm**2

       bt0inv1=1D0/bt0
       bt0inv2=1D0/bt0**2
       bt0inv3=1D0/bt0**3
       bt0inv4=1D0/bt0**4
       bt0inv5=1D0/bt0**5
       lfr=XLFR
       lqr=XLQR
       N4=5D0/3D0    !! N4 = (N^2-4)/N with N=3 for QCD
       nfv=1 !! the contribution is taken inside the flux factor 
CCCC    
c       write(*,*)"INFV=",INFV

       IF (MODL .EQ. 0 .AND. CHANNEL .EQ. 1) THEN !! DY
         ZINV = NINV2
c         FUNC  = ZINV 
       include 'resum_dy.h'
       ELSEIF(MODL .EQ. 0 .AND. CHANNEL .NE. 1) THEN
       WRITE(*,*)"Wrong channel for DY"
       STOP
c       ELSEIF(MODL.EQ.1 .OR. MODL.EQ.2)THEN !! ADD/RS
c          IF (CHANNEL .eq. 1) THEN !! qqb
c         ZINV = NINV2
c           FUNC  = ZINV 
c       include 'resum_grqqb.h'
c          ELSEIF (CHANNEL .eq. 2) THEN !! gg 
c         ZINV = NINV4
c           FUNC  = ZINV 
c       include 'resum_grgg.h'
c          ENDIF
        ENDIF

      RETURN
      END
