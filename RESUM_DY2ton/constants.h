       DOUBLE PRECISION BT0,BT1,BT2,BT3,BT4 
       DOUBLE PRECISION Bq1,Bq2,Bq3,Bq4
       DOUBLE PRECISION Dq1,Dq2,Dq3,Dq4
       DOUBLE PRECISION Dg1,Dg2,Dg3,Dg4
       DOUBLE PRECISION A1,A2,A3,A4,A5,K
       DOUBLE PRECISION Ag1,Ag2,Ag3,Ag4,Ag5
       DOUBLE PRECISION Aq1,Aq2,Aq3,Aq4,Aq5
       DOUBLE PRECISION GE,PI
       DOUBLE PRECISION Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9
c       DOUBLE PRECISION CA,CF
c       INTEGER NF

       GE=0.577215664901533D0
       PI=3.141592653589793238462643D0

       Z2=1.644934066848226D0
       Z3=1.202056903159594D0
       Z4=1.082323233711138D0
       Z5=1.036927755143370D0
       Z6=1.017343061984449D0
       Z7=1.008349277381923D0
       Z8=1.004077356197944D0
       Z9=1.002008392826082D0

!!! Beta function
       BT0 = 11./3.D0 * CA - 2./3.D0 * NF
       BT1 = 34./3.D0 * CA**2 - 10./3.D0 * CA*NF - 2.D0 * CF*NF
       BT2 = 2857./2. - 5033./18.* NF + 325./54.D0 * NF**2
       BT3 = 29243.0 - 6946.30 * NF + 405.089 * NF**2 + 1.49931 * NF**3
**** b4 is from 1701.01404 (Source file)
       BT4 = ( 8157455.0D0/16.0D0 + 621885.0D0/2.0D0*Z3
     &        - 88209.0D0/2.0D0*Z4 - 288090.0D0*Z5
     &  + NF * ( - 336460813.0D0/1944.0D0 - 4811164.0D0/81.0D0*Z3
     &           + 33935.0D0/6.0D0*Z4 + 1358995.0D0/27.0D0*Z5 )
     &  + NF**2 * ( 25960913.0D0/1944.0D0 + 698531.0D0/81.0D0*Z3
     &              - 10526.0D0/9.0D0*Z4 - 381760.0D0/81.0D0*Z5 )
     &  + NF**3 * ( - 630559.0D0/5832.0D0 - 48722.0D0/243.0D0*Z3
     &              + 1618.0D0/27.0D0*Z4 + 460.0D0/9.0D0*Z5 )
     &  + NF**4 * ( 1205.0D0/2916.0D0 - 152.0D0/81.0D0*Z3 ) )



!!! Cusp Anomalous dimensions
       A1 = 4.* CF
       K  = (67./18.D0 - Z2) * CA - 5./9.D0 * NF
       A2 = 8.* CF * K
       A3 = 16.* CF*CA**2 *(245./24.D0 - 67./9.D0 * Z2 
     1         + 11./6.D0 * Z3 + 11./5.D0 * Z2**2)
     2 + CF*CA*NF *(-836./27.D0 + 160./9.D0 * Z2 - 112./3.D0 * Z3)
     3 + CF**2*NF *(-110./3.D0 + 32.* Z3)
     4 + CF*NF**2 *(-16./27.D0)
       A4 = 20702D0 - 5171.9D0*NF + 195.5772D0*NF**2 + 3.272344D0*NF**3 !! A4 from 1805.09638 Eq.4.1
* Pade` approximation for A5
       IF(NF.eq.3)THEN
       A5 = 0.601902D0*(0.42441D0*(4.D0*PI)**5)
       ELSEIF(NF.eq.4)THEN
        A5 = 0.196796D0*(0.42441D0*(4.D0*PI)**5)
       ELSEIF(NF.eq.5)THEN
        A5 = 0.000622786D0*(0.42441D0*(4.D0*PI)**5)
       ELSE
        A5 = 0.0D0
       ENDIF

!!! Resummation Coefficients
!!! Jet function Bq for DIS 
       Bq1 = -3.* CF
       Bq2 = CF*CA  *  (-3155./54.D0 + 40.* Z3 + 44./3.D0 * Z2)
     1     + CF**2  *  (-3./2.D0 - 24.* Z3 + 12.* Z2)
     2     + CF*NF  *  (247./27.D0 - 8./3.D0 * Z2)
       Bq3 =  CF**3  * (- 29./2.D0 + 32.* Z2*Z3 - 18.* Z2
     1               - 288./5.D0 * Z2**2 - 68.* Z3 + 240.* Z5 )
     2     + CA*CF**2 * (- 46. - 16.* Z2*Z3 + 287.* Z2
     3         - 272./5.D0* Z2**2 - 712./3.D0 * Z3 - 120.* Z5 )
     4     + CA**2*CF * (- 599375./729.D0 - 176./3.D0 * Z2*Z3
     5                  + 32126./81.D0 * Z2 - 652./15.D0 * Z2**2
     6                  + 21032./27.* Z3 - 232.* Z5 )
     7     + CF**2*NF * (5501./54.D0 - 50.* Z2 + 32./9.D0 * Z3 )
     8     + CA*CF*NF* (160906./729.D0 - 9920./81.D0 * Z2
     9                  + 208./15.D0 * Z2**2 - 776./9. * Z3 )
     T     + CF*NF**2 * (- 8714./729.D0 + 232./27.D0 * Z2
     E                  - 32./27.D0 * Z3)

** Pade approximation (1/1) of B.
        IF(NF.eq.3)THEN
        Bq4 = 94704.4D0
        ELSEIF(NF.eq.4)THEN
        Bq4 = 75025.9D0
        ELSEIF(NF.eq.5)THEN
        Bq4 = 59700.0D0
        ELSE
        Bq4 = 0.0D0
        ENDIF

!!! Cusp for quarks 
      Aq1 = A1
      Aq2 = A2
      Aq3 = A3
      Aq4 = A4
      Aq5 = A5
!!!  Soft for quarks ( 0508265 )
      Dq1 = (0.D0)
      Dq2 = (CF*( CA*(-1616.D0/27.D0 + 176.D0/3.D0*Z2 + 56.D0*Z3) 
     &          + NF*(224.D0/27.D0 - 32.D0/3.D0*Z2) ) )
      Dq3 = ( CF*CA**2*(- 594058.D0/729.D0 + 98224.D0/81.D0*Z2 
     &                  + 40144.D0/27.D0*Z3 - 2992.D0/15.D0*Z2**2
     &                  - 352.D0/3.D0*Z2*Z3 - 384.D0*Z5)  
     &      + CF*CA*NF*(  125252.D0/729.D0 -29392.D0/81.D0*Z2 
     &                  - 2480.D0/9.D0*Z3 + 736.D0/15.D0*Z2**2)
     &      + CF**2*NF*(  3422.D0/27.D0 - 32.D0*Z2 
     &                  - 608.D0/9.D0*Z3 - 64.D0/5.D0*Z2**2)
     &      + CF*NF**2*(- 3712.D0/729.D0 + 640.D0/27.D0*Z2 
     &                  + 320.D0/27.D0*Z3))

!!! Cusp for gluons 
      Ag1 = CA/CF*Aq1
      Ag2 = CA/CF*Aq2
      Ag3 = CA/CF*Aq3
      Ag4 = CA/CF*Aq4
      Ag5 = CA/CF*Aq5
!!!  Soft for gluons ( 0508265 )
      Dg1 = CA/CF*Dq1
      Dg2 = CA/CF*Dq2
      Dg3 = CA/CF*Dq3

