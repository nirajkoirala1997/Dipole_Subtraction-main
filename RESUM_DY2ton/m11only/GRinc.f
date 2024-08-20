*--------------------------------------------------------------------*
      PROGRAM INCLUSIVERESUM
*--------------------------------------------------------------------*
      IMPLICIT NONE
c      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SC(11),ans(0:1)
      DATA SC/0.1D0,0.125D0,0.166667D0,0.25D0,0.5D0,1.0D0,2.0D0,4.0D0
     &       ,6.0D0,8.0D0,10.0D0/
      DATA PI/3.141592653589793238462643D0/

      INTEGER IMUR,IMUF,IMU,nkind
      INTEGER MUFI,MUFF,MURI,MURF
      INTEGER ND,NOD,NHIGGS
      INTEGER npt1,its1
      DOUBLE PRECISION answer,answer1,answer2,sd,chi2,PI,SC,ans
      INTEGER IQ,IQI,IQF,ORDI,ORDF,ISCALE,ISCALEMIN,ISCALEMAX,m11switch
      INTEGER IPDF,IPDFI,IPDFF

      DOUBLE PRECISION STEP
      DOUBLE PRECISION XMS
      DOUBLE PRECISION alphasPDF!vwgt
      DOUBLE PRECISION alphsmz,alphas
      DOUBLE PRECISION pt1


      CHARACTER fnameLO*500,fnameNLO*500,fnameNNLO*500,fnameN3LO*500
      INTEGER IMODE,NPROC
      COMMON/SVRES/IMODE

!!! Common block
      INTEGER MODL,MODE,MODP,CHANNEL,NF,ORD,INFV,counts
      DOUBLE PRECISION RS,CA,CF,TF,XNC,XMT,XMB,XMC,Q,XMUR,XMUF,AS
      DOUBLE PRECISION GF,XMZ,XMW,SW2,CW2,GMZ,GMW,ALEM,XMMH
      DOUBLE PRECISION XLFR,XLQF,XLQR
      DOUBLE PRECISION XLAMBDA,XKAP
      DOUBLE PRECISION VUD, VUS, VUB, VCD, VCS, VCB
      DOUBLE PRECISION Qmax,Qmin,Q2,XJAC,XEPS,Qbase
      CHARACTER*50 PREFIXLO,PREFIXNLO,PREFIXNNLO
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
      COMMON/VMASS/XMZ,XMW,XMMH,GMZ,GMW
      COMMON/WENBERG/SW2,CW2,ALEM
      COMMON/CKM/VUD,VUS,VUB,VCD,VCS,VCB
      COMMON/PROCESS/NPROC
      COMMON/AHIGGS/NHIGGS
      COMMON/Qvalue/Qmin,Qmax,XEPS
      COMMON/counter/counts
      COMMON/m11ONOFF/m11switch

      DOUBLE PRECISION SIGPP0(100,80,50),SIGPP1(100,80,50),
     &                 SIGPP2(100,80,50),SIGPP3(100,80,50)
      DOUBLE PRECISION error0(100,80,50),error1(100,80,50),
     &                 error2(100,80,50),error3(100,80,50)
      DOUBLE PRECISION BORN
      DOUBLE PRECISION REG
      EXTERNAL REG
      DOUBLE PRECISION BORNGRRSINC,BORNGRinc,BORNDY
      EXTERNAL BORNGRRSINC,BORNGRinc,BORNDY
      EXTERNAL alphasPDF
!!! ADD
       REAL*8 c00,ams,amh,aam1,c0,aamh
c      COMMON/CHANNEL/CHANNEL
c       COMMON/ADDRS/MODL
c       COMMON/ISET/MODE,MODP
c       common/xpar/c00,ams,amh
c       COMMON/rs_par/aam1,c0,aamh
c       COMMON/ADDPAR/XMS,NOD
****************************************************************************
      OPEN(UNIT=1,FILE="input.input",STATUS="OLD")
C     INPUT
C     -----
      READ(1,*)
      READ(1,*) RS               !  RS=ECM 
      READ(1,*)
      READ(1,*) MODP             !! MODP=1 LHC; MODP=2 TEVATRON
      READ(1,*)
      READ(1,*) MODL             ! 1=ADD, 2=RS, 3=DY 
      READ(1,*)
      READ(1,*) CHANNEL          ! 1=qqb, 2=gg 
      READ(1,*)
      READ(1,*) NPROC             ! 0=DY, 1=Z, 2=W+, 3=W-, 4=BOTH
      READ(1,*)
      READ(1,*) NHIGGS            ! Associated Higgs Production 0=OFF, 1=ON
      READ(1,*)
      READ(1,*) IMODE            ! SV=1, RESUM=2 
      READ(1,*)
      READ(1,*) ORDI, ORDF       ! order LO, NLO, NNLO 
      READ(1,*)
      READ(1,*) m11switch       ! m11 ginac 
      READ(1,*)
      READ(1,*) IQI, IQF, Qbase ,STEP       !  Q variation
      READ(1,*)
      READ(1,*) MURI, MURF       ! MuR variation
      READ(1,*)
      READ(1,*) MUFI, MUFF       ! MuF variation
      READ(1,*) 
      READ(1,*) ISCALEMIN, ISCALEMAX
      READ(1,*)
      READ(1,*)
      READ(1,*) PREFIXLO              ! LO PDF 
      READ(1,*) PREFIXNLO             ! NLO PDF 
      READ(1,*) PREFIXNNLO            ! NNLO PDF 
      READ(1,*)
      READ(1,*) IPDFI, IPDFF      ! PDFI, PDFF 
      READ(1,*)
      READ(1,*) fnameLO               ! output fine name LO
      READ(1,*) fnameNLO              ! output fine name NLO
      READ(1,*) fnameNNLO             ! output fine name NNLO
      READ(1,*) fnameN3LO             ! output fine name N3LO
      READ(1,*)
      CLOSE(1)
**************************INFO**************************************************
!!! Vegas info
      OPEN(unit=10,file='vegas.input',status='unknown')
      READ (10,*) pt1          ! vegas points     LO 2 body
      npt1 = pt1
      READ (10,*) its1          ! vegas iterations LO 2 body
      CLOSE(10)
      ND = 5
****************************************************************************
c        OPEN(2,file='in_add.input')
c        READ(2,*) XMS          !M_S Fundamental Planck scale
c        READ(2,*) NOD          !d:[2,6] extra dim
c        CLOSE(2)
c        OPEN(3,file='in_rs.input')
c        READ(3,*) aam1       !M1 1st excited RS-KK mode
c        READ(3,*) C0         !c0 effective RS coupling
c        READ(3,*) aamh       !XMH
c        CLOSE(3)
****************************************************************************

      GF   = 1.16637900E-05
      XMZ  = 91.1876D0
      XMW  = 80.3850d0
      GMZ  = 2.4952d0
      GMW  = 2.0850d0
      CW2  = XMW**2.0d0/XMZ**2.0d0
      SW2  = 1.0d0 - CW2
c      ALEM = GF * (8.0d0*SW2*CW2*XMZ**2.0d0)/(4.d0*sqrt(2.d0)*PI)
      ALEM = 1.0d0/128.0
      XMT  = 172.760D0
      XMMH  = 125.10d0

c-------------------------------------------------------------------

       VUD = 0.97446D0
       VUS = 0.22452D0
       VUB = 0.00365D0
       VCD = 0.22438D0
       VCS = 0.97359D0
       VCB = 0.04214D0

      write(*,*)'******************************************************'
      IF (NPROC.EQ.0) THEN
      WRITE(*,*)'INVARIANT MASS DISTRIBUTION OF ZZ AT THE LHC '
      write(*,*)'******************************************************'
      ENDIF
      WRITE(*,*)'CENTRE OF ENERGY =',RS, 'TEV'
      WRITE(*,*)'MASS OF HIGGS =',XMMH, 'MASS OF Z =',XMZ
      WRITE(*,*)'FERMI CONSTANT =',GF, 'FINE STRUCTURE CONSTANT=',1.0D0/
     . ALEM
      WRITE(*,*)'SIN^2(THETA_w) =',SW2, 'COS^2(THETA_w)',CW2
      write(*,*)'******************************************************'
      
c*******************************************
      ! ORDER VARIATION !
c*******************************************
      DO ORD = ORDI, ORDF  


c **************************************
!!!     PDF from LHAPDF 
c **************************************
C     -----------------------------
c      WRITE(*,*)"ORD=",ORD

      IF     (ORD.EQ.0) THEN
      PREFIX = PREFIXLO
      ELSEIF (ORD.EQ.1) THEN
      PREFIX = PREFIXNLO
      ELSEIF (ORD.EQ.2) THEN
      PREFIX = PREFIXNNLO
      ELSEIF (ORD.EQ.3) THEN
      PREFIX = PREFIXNNLO
      ENDIF

      call InitPDFsetByName(PREFIX)

c*******************************************
      ! PDF VARIATION !
c*******************************************

      DO IPDF = IPDFI, IPDFF   !! PDF Variation
      call initpdf(IPDF)

!!! Colors
      CA  = 3.0D0
      CF  = 4.0D0/3.0D0
      TF  = 1.0D0/2.0D0
      XNC = 3.0D0
      NF  = 5.0D0
c      Qbase = 175d0 

      DO IQ=IQI, IQF  !! Q variation
      Q =  Qbase + STEP*IQ    


c      Q = 500.0d0
      XEPS = 5d0
      Qmin = Q - XEPS
      Qmax = Q + XEPS

       write(*,*)'Qmin, Qmax = ',Qmin, Qmax
c      DO IMUR=MURI,MURF       !! MuR Variation
c      DO IMUF=MUFI,MUFF       !! MuF Variation
      DO ISCALE = ISCALEMIN,ISCALEMAX
      IF (ISCALE .EQ. 1) THEN
        XMUR = 1.0D0/2.0D0*Q
        XMUF = 1.0D0/2.0D0*Q
c        WRITE(*,*)'***********(1/2,1/2)**********'
      ELSEIF (ISCALE .EQ. 2) THEN
        XMUR = 1.0D0/2.0D0*Q
        XMUF = Q
c        WRITE(*,*)'***********(1/2,1)**********'
      ELSEIF (ISCALE .EQ. 3) THEN
        XMUR = Q
        XMUF = 1.0D0/2.0D0*Q
c        WRITE(*,*)'***********(1,1/2)**********'
      ELSEIF (ISCALE .EQ. 4) THEN
        XMUR = Q
        XMUF = Q
c        WRITE(*,*)'***********(1,1)**********'
      ELSEIF (ISCALE .EQ. 5) THEN
        XMUR = Q
        XMUF = 2.0D0*Q
c        WRITE(*,*)'***********(1,2)**********'
      ELSEIF (ISCALE .EQ. 6) THEN
        XMUR = 2.0d0*Q
        XMUF = Q
c        WRITE(*,*)'***********(2,1)**********'
      ELSEIF (ISCALE .EQ. 7) THEN
        XMUR = 2.0d0*Q
        XMUF = 2.0d0*Q
c        WRITE(*,*)'***********(2,2)**********'
      ENDIF

       
!!! AlphaS/4/Pi
C      AS = alphasPDF(XMUR)/4.D0/PI
      alphsmz = alphasPDF(XMZ)
      nkind = 2    !! what is this ?
      AS = alphas(nkind,ord,xmz,prefix,alphsmz,xmur)/4.0d0/PI


      XLFR = DLOG(XMUF**2/XMUR**2)
      XLQF = DLOG(Q**2/XMUF**2)
      XLQR = DLOG(Q**2/XMUR**2)
c      WRITE(*,*)">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*)'******************************************************'
      WRITE(*,*)'Running Order',ORD,'with PDF=',PREFIX
      WRITE(*,*)"Q=",Q,"MUR=",XMUR,"MUF=",XMUF
      write(*,*)'******************************************************'


CCC------------------------------VEGAS-------------------------------CCC
       IF (NPROC.EQ.0 .OR. NPROC .EQ. 1) THEN
       IF (MODL.EQ.0 .AND. ORD.EQ.3)THEN
        DO INFV=0,1
        counts = 0
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(ND,npt1,its1,REG,answer1,sd,chi2) !! all other cases
        ans(INFV) = answer1
        ENDDO
        answer=ans(0)+ans(1)
       ELSE
        counts = 0
        call brm48i(40,0,0) ! initialize random number generator
        call vsup(ND,npt1,its1,REG,answer,sd,chi2)
       ENDIF
       ENDIF
       write(*,*)"answer , counts=",answer , counts
       write(*,*)"error =",sd

!!! 
C------------------------------- LO--------------------------------------------C
      IF     (ORD.EQ.0) THEN
        SIGPP0(IPDF,IQ,ISCALE)=answer
        error0(IPDF,IQ,ISCALE)=sd
        WRITE(*,*)'XSECT-00 =',Q,XMUR,XMUF,SIGPP0(IPDF,IQ,ISCALE),
     &                                  error0(IPDF,IQ,ISCALE)
C--------------------- NLO--------------------------------------------C
      ELSEIF (ORD.EQ.1) THEN
        SIGPP1(IPDF,IQ,ISCALE)=answer
        error1(IPDF,IQ,ISCALE)=sd
        WRITE(*,*)'XSECT-01 =',Q,XMUR,XMUF,SIGPP1(IPDF,IQ,ISCALE),
     &                                  error1(IPDF,IQ,ISCALE)
C------------------------------- NNLO--------------------------------------------C
      ELSEIF (ORD.EQ.2) THEN
        SIGPP2(IPDF,IQ,ISCALE)=answer
        error2(IPDF,IQ,ISCALE)=sd
        WRITE(*,*)'XSECT-02 =',Q,XMUR,XMUF,SIGPP2(IPDF,IQ,ISCALE),
     &                                  error2(IPDF,IQ,ISCALE)
C------------------ NNNLO--------------------------------------------C
      ELSEIF (ORD.EQ.3) THEN
        SIGPP3(IPDF,IQ,ISCALE)=answer
        error3(IPDF,IQ,ISCALE)=sd
        WRITE(*,*)'XSECT-03 =',Q,XMUR,XMUF,SIGPP3(IPDF,IQ,ISCALE),
     &                                  error3(IPDF,IQ,ISCALE)
      ENDIF

       ENDDO !! 7-point
       ENDDO !! Q   LOOP
       ENDDO !! PDF LOOP
       ENDDO !! ORD LOOP
            
!!!!---------------- writing into files (same structure) ----------------------!!!!
      DO ORD = ORDI, ORDF      !! ORDER Variation
      DO IPDF = IPDFI, IPDFF   !! PDF Variation
 
       IF     (ORD.EQ.0) THEN
        IF (IPDF .NE. 0) THEN !! meaning pdf variation
        write(fnameLO,'("lores_",i3.3,".dat")') iPDF !! for PDF variation
        ENDIF
        OPEN (unit = iPDF, file = fnameLO)
       ELSEIF (ORD.EQ.1) THEN
        IF (IPDF .NE. 0) THEN !! meaning pdf variation
        write(fnameNLO,'("nlores_",i3.3,".dat")') iPDF !! for PDF variation
        ENDIF
        OPEN (unit = iPDF, file = fnameNLO)
       ELSEIF (ORD.EQ.2) THEN
        IF (IPDF .NE. 0) THEN !! meaning pdf variation
        write(fnameNNLO,'("nnlores_",i3.3,".dat")') iPDF !! for PDF variation
        ENDIF
        OPEN (unit = iPDF, file = fnameNNLO)
       ELSEIF (ORD.EQ.3) THEN
        IF (IPDF .NE. 0) THEN !! meaning pdf variation
        write(fnameN3LO,'("n3lores_",i3.3,".dat")') iPDF !! for PDF variation
        ENDIF
        OPEN (unit = iPDF, file = fnameN3LO)
       ENDIF


c      Qbase = 175d0 

      DO IQ=IQI, IQF  !! Q variation

      Q =  Qbase + STEP*IQ    

c      Q =  STEP*IQ
      DO ISCALE = ISCALEMIN,ISCALEMAX      !! 7-point scale variation
      IF (ISCALE .EQ. 1) THEN
        XMUR = 1.0D0/2.0D0*Q
        XMUF = 1.0D0/2.0D0*Q
c        WRITE(*,*)'***********(1/2,1/2)**********'
      ELSEIF (ISCALE .EQ. 2) THEN
        XMUR = 1.0D0/2.0D0*Q
        XMUF = Q
c        WRITE(*,*)'***********(1/2,1)**********'
      ELSEIF (ISCALE .EQ. 3) THEN
        XMUR = Q
        XMUF = 1.0D0/2.0D0*Q
c        WRITE(*,*)'***********(1,1/2)**********'
      ELSEIF (ISCALE .EQ. 4) THEN
        XMUR = Q
        XMUF = Q
c        WRITE(*,*)'***********(1,1)**********'
      ELSEIF (ISCALE .EQ. 5) THEN
        XMUR = Q
        XMUF = 2.0D0*Q
c        WRITE(*,*)'***********(1,2)**********'
      ELSEIF (ISCALE .EQ. 6) THEN
        XMUR = 2.0d0*Q
        XMUF = Q
c        WRITE(*,*)'***********(2,1)**********'
      ELSEIF (ISCALE .EQ. 7) THEN
        XMUR = 2.0d0*Q
        XMUF = 2.0d0*Q
c        WRITE(*,*)'***********(2,2)**********'
      ENDIF
      IF     (ORD.EQ.0) THEN
        WRITE(iPDF,*)Q,XMUR,XMUF,SIGPP0(IPDF,IQ,ISCALE),
     &                          error0(IPDF,IQ,ISCALE) 
      ELSEIF (ORD.EQ.1) THEN
        WRITE(iPDF,*)Q,XMUR,XMUF,SIGPP1(IPDF,IQ,ISCALE),
     &                          error1(IPDF,IQ,ISCALE) 
      ELSEIF (ORD.EQ.2) THEN
        WRITE(iPDF,*)Q,XMUR,XMUF,SIGPP2(IPDF,IQ,ISCALE),
     &                          error2(IPDF,IQ,ISCALE) 
      ELSEIF (ORD.EQ.3) THEN
        WRITE(iPDF,*)Q,XMUR,XMUF,SIGPP3(IPDF,IQ,ISCALE),
     &                          error3(IPDF,IQ,ISCALE) 
      ENDIF

      ENDDO !! 7-point
      ENDDO !! Q   LOOP
       CLOSE(unit=iPDF)
      ENDDO !! PDF LOOP
      ENDDO !! ORD LOOP

       RETURN
       END
