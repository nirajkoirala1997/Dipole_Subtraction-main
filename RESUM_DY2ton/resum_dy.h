***** (**w = 2*bt0*a(muR^2)*Log[Nb]=2*lam**)  ***** 
*****   Log1mw = Log[1-w] ***** 

        IF (ORD .EQ. 0) THEN
*****   LL ***** 

       gll=  + lnNb * ( 8*bt0inv1*Cf )
      gll = gll + lnNb*zlwm * (  - 8*bt0inv1*Cf )
      gll = gll + lnNb*winv*zlwm * ( 8*bt0inv1*Cf ) 

       CHg0tll=  + 1 

       foll=  + 1 

CC>>>>>>>>>>>>>>>>>>>>>>>>>LO
      IF(IMODE .EQ. 1) THEN
       FUNC = foll*Zinv                   !!SV
      ELSEIF(IMODE .EQ. 2) THEN
         IF (REAL(w).ge. 1.0D0) THEN
         FUNC =(0.0D0,0.0D0)
         ELSE
         FUNC = (CHg0tll*zexp(gll)-foll)*Zinv   !!RESUM-SV
         ENDIF
      ENDIF

        ELSEIF (ORD .EQ. 1) THEN
*****   NLL ***** 

       gnll=  + zlwm2 * (  - 4*bt0inv3*Cf**2*nf - 20.D0/3.D0*bt0inv3*Ca
     &    *Cf*nf + 68.D0/3.D0*bt0inv3*Ca**2*Cf )
      gnll = gnll + zlwm * (  - 8*bt0inv3*Cf**2*nf - 40.D0/3.D0*bt0inv3
     &    *Ca*Cf*nf + 136.D0/3.D0*bt0inv3*Ca**2*Cf + 40.D0/9.D0*bt0inv2
     &    *Cf*nf - 268.D0/9.D0*bt0inv2*Ca*Cf + 8*z2*bt0inv2*Ca*Cf )
      gnll = gnll + w * (  - 8*bt0inv3*Cf**2*nf - 40.D0/3.D0*bt0inv3*Ca
     &    *Cf*nf + 136.D0/3.D0*bt0inv3*Ca**2*Cf + 40.D0/9.D0*bt0inv2*Cf
     &    *nf - 268.D0/9.D0*bt0inv2*Ca*Cf + 8*z2*bt0inv2*Ca*Cf )
      gnll = gnll + lqr*zlwm * ( 4*bt0inv1*Cf )
      gnll = gnll + lfr*w * ( 4*bt0inv1*Cf )
      gnll = gnll + lnNb * ( 8*bt0inv1*Cf )
      gnll = gnll + lnNb*zlwm * (  - 8*bt0inv1*Cf )
      gnll = gnll + lnNb*winv*zlwm * ( 8*bt0inv1*Cf ) 

c      CHg0tnll=  + as * (  - 16*Cf + 16*z2*Cf )
c      CHg0tnll = CHg0tnll + lqr*as * ( 6*Cf )
c      CHg0tnll = CHg0tnll + lfr*as * (  - 6*Cf )
      

      CHg0tnll =  virt1
      CHg0tnll = CHg0tnll + 1 
      
       fonll = 0.0d0
c      fonll = virt1
c      fonll = fonll + 1 

c       fonll = fonll + as * (2.d0 * Cf * z2)
c       fonll = fonll + lqr*as * ( 6*Cf )
c       fonll = fonll + lfr*as * (  - 6*Cf )
c       fonll = fonll + lnNb*lqr*as * (  - 8*Cf )
c       fonll = fonll + lnNb*lfr*as * ( 8*Cf )
       fonll = fonll + lnNb**2*as * ( 8*Cf )

 
c      fonll = fonll + lqr*as * ( 6*Cf )
c      fonll = fonll + lfr*as * (  - 6*Cf )
c      fonll = fonll + lnNb*lqr*as * (  - 8*Cf )
c      fonll = fonll + lnNb*lfr*as * ( 8*Cf )
c      fonll = fonll + lnNb**2*as * ( 8*Cf )
c
c      fonll=  + as * (  - 16*Cf + 16*z2*Cf )
c      fonll = fonll + lqr*as * ( 6*Cf )
c      fonll = fonll + lfr*as * (  - 6*Cf )
cc      fonll = fonll + lnNb*lqr*as * (  - 8*Cf )
c      fonll = fonll + lnNb*lfr*as * ( 8*Cf )
c      fonll = fonll + lnNb**2*as * ( 8*Cf )
      fonll = fonll + 1 

CC>>>>>>>>>>>>>>>>>>>>>>>>>NLO
      IF(IMODE .EQ. 1) THEN
       FUNC = fonll*Zinv                   !!SV
      ELSEIF(IMODE .EQ. 2) THEN
         IF (REAL(w).ge. 1.0D0) THEN
         FUNC =(0.0D0,0.0D0)
         ELSE
         FUNC = (CHg0tnll*zexp(gnll)-fonll)*Zinv   !!RESUM-SV
         ENDIF
      ENDIF

        ELSEIF (ORD .EQ. 2) THEN
*****   NNLL ***** 

       gnnll=  + zlwm2 * (  - 4*bt0inv3*Cf**2*nf - 20.D0/3.D0*bt0inv3*
     &    Ca*Cf*nf + 68.D0/3.D0*bt0inv3*Ca**2*Cf )
      gnnll = gnnll + zlwm * (  - 8*bt0inv3*Cf**2*nf - 40.D0/3.D0*
     &    bt0inv3*Ca*Cf*nf + 136.D0/3.D0*bt0inv3*Ca**2*Cf + 40.D0/9.D0*
     &    bt0inv2*Cf*nf - 268.D0/9.D0*bt0inv2*Ca*Cf + 8*z2*bt0inv2*Ca*
     &    Cf )
      gnnll = gnnll + w * (  - 8*bt0inv3*Cf**2*nf - 40.D0/3.D0*bt0inv3*
     &    Ca*Cf*nf + 136.D0/3.D0*bt0inv3*Ca**2*Cf + 40.D0/9.D0*bt0inv2*
     &    Cf*nf - 268.D0/9.D0*bt0inv2*Ca*Cf + 8*z2*bt0inv2*Ca*Cf )
      gnnll = gnnll + as*wminv*zlwm2 * ( 8*bt0inv4*Cf**3*nf**2 + 80.D0/
     &    3.D0*bt0inv4*Ca*Cf**2*nf**2 + 200.D0/9.D0*bt0inv4*Ca**2*Cf*
     &    nf**2 - 272.D0/3.D0*bt0inv4*Ca**2*Cf**2*nf - 1360.D0/9.D0*
     &    bt0inv4*Ca**3*Cf*nf + 2312.D0/9.D0*bt0inv4*Ca**4*Cf )
      gnnll = gnnll + as*wminv*zlwm * (  - 4*bt0inv3*Cf**2*nf**2 + 4*
     &    bt0inv3*Cf**3*nf - 242.D0/27.D0*bt0inv3*Ca*Cf*nf**2 + 14*
     &    bt0inv3*Ca*Cf**2*nf + 1210.D0/27.D0*bt0inv3*Ca**2*Cf*nf - 
     &    3398.D0/27.D0*bt0inv3*Ca**3*Cf - 16*z2*bt0inv3*Ca*Cf**2*nf - 
     &    80.D0/3.D0*z2*bt0inv3*Ca**2*Cf*nf + 272.D0/3.D0*z2*bt0inv3*
     &    Ca**3*Cf )
      gnnll = gnnll + as*w*wminv * (  - 4*bt0inv3*Cf**2*nf**2 + 4*
     &    bt0inv3*Cf**3*nf - 242.D0/27.D0*bt0inv3*Ca*Cf*nf**2 + 14*
     &    bt0inv3*Ca*Cf**2*nf + 1210.D0/27.D0*bt0inv3*Ca**2*Cf*nf - 
     &    3398.D0/27.D0*bt0inv3*Ca**3*Cf - 112.D0/27.D0*bt0inv1*Cf*nf
     &     + 808.D0/27.D0*bt0inv1*Ca*Cf - 28*z3*bt0inv1*Ca*Cf + 8*z2*Cf
     &     - 16*z2*bt0inv3*Ca*Cf**2*nf - 80.D0/3.D0*z2*bt0inv3*Ca**2*Cf
     &    *nf + 272.D0/3.D0*z2*bt0inv3*Ca**3*Cf + 16.D0/3.D0*z2*bt0inv1
     &    *Cf*nf - 88.D0/3.D0*z2*bt0inv1*Ca*Cf )
      gnnll = gnnll + as*w*wminv*zlwm * ( 16*bt0inv4*Cf**3*nf**2 + 160.D
     &    0/3.D0*bt0inv4*Ca*Cf**2*nf**2 + 400.D0/9.D0*bt0inv4*Ca**2*Cf*
     &    nf**2 - 544.D0/3.D0*bt0inv4*Ca**2*Cf**2*nf - 2720.D0/9.D0*
     &    bt0inv4*Ca**3*Cf*nf + 4624.D0/9.D0*bt0inv4*Ca**4*Cf - 44.D0/9.
     &    D0*bt0inv3*Cf**2*nf**2 - 4*bt0inv3*Cf**3*nf - 158.D0/27.D0*
     &    bt0inv3*Ca*Cf*nf**2 + 410.D0/9.D0*bt0inv3*Ca*Cf**2*nf + 2830.D
     &    0/27.D0*bt0inv3*Ca**2*Cf*nf - 5714.D0/27.D0*bt0inv3*Ca**3*Cf
     &     )
      gnnll = gnnll + as*w**2*wminv * ( 8*bt0inv4*Cf**3*nf**2 + 80.D0/3.
     &    D0*bt0inv4*Ca*Cf**2*nf**2 + 200.D0/9.D0*bt0inv4*Ca**2*Cf*
     &    nf**2 - 272.D0/3.D0*bt0inv4*Ca**2*Cf**2*nf - 1360.D0/9.D0*
     &    bt0inv4*Ca**3*Cf*nf + 2312.D0/9.D0*bt0inv4*Ca**4*Cf - 62.D0/9.
     &    D0*bt0inv3*Cf**2*nf**2 - 2*bt0inv3*Cf**3*nf - 31.D0/3.D0*
     &    bt0inv3*Ca*Cf*nf**2 + 473.D0/9.D0*bt0inv3*Ca*Cf**2*nf + 1145.D
     &    0/9.D0*bt0inv3*Ca**2*Cf*nf - 2471.D0/9.D0*bt0inv3*Ca**3*Cf - 
     &    8.D0/27.D0*bt0inv2*Cf*nf**2 - 55.D0/3.D0*bt0inv2*Cf**2*nf - 
     &    418.D0/27.D0*bt0inv2*Ca*Cf*nf + 245.D0/3.D0*bt0inv2*Ca**2*Cf
     &     + 16*z3*bt0inv2*Cf**2*nf - 56.D0/3.D0*z3*bt0inv2*Ca*Cf*nf + 
     &    44.D0/3.D0*z3*bt0inv2*Ca**2*Cf - 8*z2*bt0inv3*Ca*Cf**2*nf - 
     &    40.D0/3.D0*z2*bt0inv3*Ca**2*Cf*nf + 136.D0/3.D0*z2*bt0inv3*
     &    Ca**3*Cf + 80.D0/9.D0*z2*bt0inv2*Ca*Cf*nf - 536.D0/9.D0*z2*
     &    bt0inv2*Ca**2*Cf + 88.D0/5.D0*z2**2*bt0inv2*Ca**2*Cf )
      gnnll = gnnll + lqr*zlwm * ( 4*bt0inv1*Cf )
      gnnll = gnnll + lqr*as*wminv*zlwm * (  - 8*bt0inv2*Cf**2*nf - 40.D
     &    0/3.D0*bt0inv2*Ca*Cf*nf + 136.D0/3.D0*bt0inv2*Ca**2*Cf )
      gnnll = gnnll + lqr*as*w*wminv * (  - 8*bt0inv2*Cf**2*nf - 40.D0/
     &    3.D0*bt0inv2*Ca*Cf*nf + 136.D0/3.D0*bt0inv2*Ca**2*Cf + 40.D0/
     &    9.D0*bt0inv1*Cf*nf - 268.D0/9.D0*bt0inv1*Ca*Cf + 8*z2*bt0inv1
     &    *Ca*Cf )
      gnnll = gnnll + lqr**2*as*w*wminv * ( 2*Cf )
      gnnll = gnnll + lfr*w * ( 4*bt0inv1*Cf )
      gnnll = gnnll + lfr*as*w*wminv * (  - 40.D0/9.D0*bt0inv1*Cf*nf + 
     &    268.D0/9.D0*bt0inv1*Ca*Cf - 8*z2*bt0inv1*Ca*Cf )
      gnnll = gnnll + lfr*as*w**2*wminv * ( 40.D0/9.D0*bt0inv1*Cf*nf - 
     &    268.D0/9.D0*bt0inv1*Ca*Cf + 8*z2*bt0inv1*Ca*Cf )
      gnnll = gnnll + lfr**2*as*w*wminv * (  - 2*Cf )
      gnnll = gnnll + lfr**2*as*w**2*wminv * ( 2*Cf )
      gnnll = gnnll + lnNb * ( 8*bt0inv1*Cf )
      gnnll = gnnll + lnNb*zlwm * (  - 8*bt0inv1*Cf )
      gnnll = gnnll + lnNb*winv*zlwm * ( 8*bt0inv1*Cf ) 

      CHg0tnnll = virt1
      CHg0tnnll = CHg0tnnll + 1 

      fonnll = virt1
      fonnll = fonnll + 1 
CC>>>>>>>>>>>>>>>>>>>>>>>>>NNLO
      IF(IMODE .EQ. 1) THEN
       FUNC = fonnll*Zinv                   !!SV
      ELSEIF(IMODE .EQ. 2) THEN
         IF (REAL(w).ge. 1.0D0) THEN
         FUNC =(0.0D0,0.0D0)
         ELSE
         FUNC = (CHg0tnnll*zexp(gnnll)-fonnll)*Zinv   !!RESUM-SV
         ENDIF
      ENDIF

        ELSEIF (ORD .EQ. 3) THEN
*****   NNNLL ***** 

       gnnnll=  + zlwm2 * (  - 4*bt0inv3*Cf**2*nf - 20.D0/3.D0*bt0inv3*
     &    Ca*Cf*nf + 68.D0/3.D0*bt0inv3*Ca**2*Cf )
      gnnnll = gnnnll + zlwm * (  - 8*bt0inv3*Cf**2*nf - 40.D0/3.D0*
     &    bt0inv3*Ca*Cf*nf + 136.D0/3.D0*bt0inv3*Ca**2*Cf + 40.D0/9.D0*
     &    bt0inv2*Cf*nf - 268.D0/9.D0*bt0inv2*Ca*Cf + 8*z2*bt0inv2*Ca*
     &    Cf )
      gnnnll = gnnnll + w * (  - 8*bt0inv3*Cf**2*nf - 40.D0/3.D0*
     &    bt0inv3*Ca*Cf*nf + 136.D0/3.D0*bt0inv3*Ca**2*Cf + 40.D0/9.D0*
     &    bt0inv2*Cf*nf - 268.D0/9.D0*bt0inv2*Ca*Cf + 8*z2*bt0inv2*Ca*
     &    Cf )
      gnnnll = gnnnll + as*wminv*zlwm2 * ( 8*bt0inv4*Cf**3*nf**2 + 80.D0
     &    /3.D0*bt0inv4*Ca*Cf**2*nf**2 + 200.D0/9.D0*bt0inv4*Ca**2*Cf*
     &    nf**2 - 272.D0/3.D0*bt0inv4*Ca**2*Cf**2*nf - 1360.D0/9.D0*
     &    bt0inv4*Ca**3*Cf*nf + 2312.D0/9.D0*bt0inv4*Ca**4*Cf )
      gnnnll = gnnnll + as*wminv*zlwm * (  - 4*bt0inv3*Cf**2*nf**2 + 4*
     &    bt0inv3*Cf**3*nf - 242.D0/27.D0*bt0inv3*Ca*Cf*nf**2 + 14*
     &    bt0inv3*Ca*Cf**2*nf + 1210.D0/27.D0*bt0inv3*Ca**2*Cf*nf - 
     &    3398.D0/27.D0*bt0inv3*Ca**3*Cf - 16*z2*bt0inv3*Ca*Cf**2*nf - 
     &    80.D0/3.D0*z2*bt0inv3*Ca**2*Cf*nf + 272.D0/3.D0*z2*bt0inv3*
     &    Ca**3*Cf )
      gnnnll = gnnnll + as*w*wminv * (  - 4*bt0inv3*Cf**2*nf**2 + 4*
     &    bt0inv3*Cf**3*nf - 242.D0/27.D0*bt0inv3*Ca*Cf*nf**2 + 14*
     &    bt0inv3*Ca*Cf**2*nf + 1210.D0/27.D0*bt0inv3*Ca**2*Cf*nf - 
     &    3398.D0/27.D0*bt0inv3*Ca**3*Cf - 112.D0/27.D0*bt0inv1*Cf*nf
     &     + 808.D0/27.D0*bt0inv1*Ca*Cf - 28*z3*bt0inv1*Ca*Cf + 8*z2*Cf
     &     - 16*z2*bt0inv3*Ca*Cf**2*nf - 80.D0/3.D0*z2*bt0inv3*Ca**2*Cf
     &    *nf + 272.D0/3.D0*z2*bt0inv3*Ca**3*Cf + 16.D0/3.D0*z2*bt0inv1
     &    *Cf*nf - 88.D0/3.D0*z2*bt0inv1*Ca*Cf )
      gnnnll = gnnnll + as*w*wminv*zlwm * ( 16*bt0inv4*Cf**3*nf**2 + 
     &    160.D0/3.D0*bt0inv4*Ca*Cf**2*nf**2 + 400.D0/9.D0*bt0inv4*
     &    Ca**2*Cf*nf**2 - 544.D0/3.D0*bt0inv4*Ca**2*Cf**2*nf - 2720.D0/
     &    9.D0*bt0inv4*Ca**3*Cf*nf + 4624.D0/9.D0*bt0inv4*Ca**4*Cf - 44.
     &    D0/9.D0*bt0inv3*Cf**2*nf**2 - 4*bt0inv3*Cf**3*nf - 158.D0/27.D
     &    0*bt0inv3*Ca*Cf*nf**2 + 410.D0/9.D0*bt0inv3*Ca*Cf**2*nf + 
     &    2830.D0/27.D0*bt0inv3*Ca**2*Cf*nf - 5714.D0/27.D0*bt0inv3*
     &    Ca**3*Cf )
      gnnnll = gnnnll + as*w**2*wminv * ( 8*bt0inv4*Cf**3*nf**2 + 80.D0/
     &    3.D0*bt0inv4*Ca*Cf**2*nf**2 + 200.D0/9.D0*bt0inv4*Ca**2*Cf*
     &    nf**2 - 272.D0/3.D0*bt0inv4*Ca**2*Cf**2*nf - 1360.D0/9.D0*
     &    bt0inv4*Ca**3*Cf*nf + 2312.D0/9.D0*bt0inv4*Ca**4*Cf - 62.D0/9.
     &    D0*bt0inv3*Cf**2*nf**2 - 2*bt0inv3*Cf**3*nf - 31.D0/3.D0*
     &    bt0inv3*Ca*Cf*nf**2 + 473.D0/9.D0*bt0inv3*Ca*Cf**2*nf + 1145.D
     &    0/9.D0*bt0inv3*Ca**2*Cf*nf - 2471.D0/9.D0*bt0inv3*Ca**3*Cf - 
     &    8.D0/27.D0*bt0inv2*Cf*nf**2 - 55.D0/3.D0*bt0inv2*Cf**2*nf - 
     &    418.D0/27.D0*bt0inv2*Ca*Cf*nf + 245.D0/3.D0*bt0inv2*Ca**2*Cf
     &     + 16*z3*bt0inv2*Cf**2*nf - 56.D0/3.D0*z3*bt0inv2*Ca*Cf*nf + 
     &    44.D0/3.D0*z3*bt0inv2*Ca**2*Cf - 8*z2*bt0inv3*Ca*Cf**2*nf - 
     &    40.D0/3.D0*z2*bt0inv3*Ca**2*Cf*nf + 136.D0/3.D0*z2*bt0inv3*
     &    Ca**3*Cf + 80.D0/9.D0*z2*bt0inv2*Ca*Cf*nf - 536.D0/9.D0*z2*
     &    bt0inv2*Ca**2*Cf + 88.D0/5.D0*z2**2*bt0inv2*Ca**2*Cf )
      gnnnll = gnnnll + as**2*wminv2*zlwm2 * (  - 80.D0/9.D0*bt0inv4*
     &    Cf**3*nf**3 - 800.D0/27.D0*bt0inv4*Ca*Cf**2*nf**3 + 536.D0/9.D
     &    0*bt0inv4*Ca*Cf**3*nf**2 - 2000.D0/81.D0*bt0inv4*Ca**2*Cf*
     &    nf**3 + 8080.D0/27.D0*bt0inv4*Ca**2*Cf**2*nf**2 + 1000.D0/3.D0
     &    *bt0inv4*Ca**3*Cf*nf**2 - 18224.D0/27.D0*bt0inv4*Ca**3*Cf**2*
     &    nf - 38080.D0/27.D0*bt0inv4*Ca**4*Cf*nf + 154904.D0/81.D0*
     &    bt0inv4*Ca**5*Cf - 16*z2*bt0inv4*Ca*Cf**3*nf**2 - 160.D0/3.D0
     &    *z2*bt0inv4*Ca**2*Cf**2*nf**2 - 400.D0/9.D0*z2*bt0inv4*Ca**3*
     &    Cf*nf**2 + 544.D0/3.D0*z2*bt0inv4*Ca**3*Cf**2*nf + 2720.D0/9.D
     &    0*z2*bt0inv4*Ca**4*Cf*nf - 4624.D0/9.D0*z2*bt0inv4*Ca**5*Cf )
      gnnnll = gnnnll + as**2*wminv2*zlwm * (  - 88.D0/3.D0*CA**(-2)*
     &    bt0inv3*Cf*nf**2 + 64*CA**(-2)*z3*bt0inv3*Cf*nf**2 - 4*
     &    bt0inv4*Cf**3*nf**3 + 4*bt0inv4*Cf**4*nf**2 - 422.D0/27.D0*
     &    bt0inv4*Ca*Cf**2*nf**3 + 62.D0/3.D0*bt0inv4*Ca*Cf**3*nf**2 - 
     &    1210.D0/81.D0*bt0inv4*Ca**2*Cf*nf**3 + 2452.D0/27.D0*bt0inv4*
     &    Ca**2*Cf**2*nf**2 - 68.D0/3.D0*bt0inv4*Ca**2*Cf**3*nf + 3388.D
     &    0/27.D0*bt0inv4*Ca**3*Cf*nf**2 - 5540.D0/27.D0*bt0inv4*Ca**3*
     &    Cf**2*nf - 12520.D0/27.D0*bt0inv4*Ca**4*Cf*nf + 57766.D0/81.D0
     &    *bt0inv4*Ca**5*Cf + 88.D0/9.D0*bt0inv3*Cf*nf**2 + 164.D0/243.D
     &    0*bt0inv3*Cf**2*nf**3 - 314.D0/27.D0*bt0inv3*Cf**3*nf**2 + 46
     &    *bt0inv3*Cf**4*nf + 128.D0/9.D0*bt0inv3*Ca*Cf*nf - 134.D0/243.
     &    D0*bt0inv3*Ca*Cf*nf**3 - 13798.D0/243.D0*bt0inv3*Ca*Cf**2*
     &    nf**2 - 4204.D0/27.D0*bt0inv3*Ca*Cf**3*nf - 80.D0/3.D0*
     &    bt0inv3*Ca**2*Cf - 25.D0/27.D0*bt0inv3*Ca**2*Cf*nf**2 + 97253.
     &    D0/243.D0*bt0inv3*Ca**2*Cf**2*nf - 2689.D0/81.D0*bt0inv3*
     &    Ca**3*Cf*nf )
      gnnnll = gnnnll + as**2*wminv2*zlwm * (  - 74437.D0/243.D0*
     &    bt0inv3*Ca**4*Cf - 224.D0/27.D0*bt0inv2*Cf**2*nf**2 - 1120.D0/
     &    81.D0*bt0inv2*Ca*Cf*nf**2 + 1616.D0/27.D0*bt0inv2*Ca*Cf**2*nf
     &     + 11888.D0/81.D0*bt0inv2*Ca**2*Cf*nf - 27472.D0/81.D0*
     &    bt0inv2*Ca**3*Cf - 64.D0/3.D0*z3*bt0inv3*Cf*nf**2 - 64.D0/9.D0
     &    *z3*bt0inv3*Cf**3*nf**2 - 416.D0/3.D0*z3*bt0inv3*Ca*Cf*nf + 
     &    368.D0/9.D0*z3*bt0inv3*Ca*Cf**2*nf**2 + 352.D0/9.D0*z3*
     &    bt0inv3*Ca*Cf**3*nf + 704*z3*bt0inv3*Ca**2*Cf - 416.D0/9.D0*
     &    z3*bt0inv3*Ca**2*Cf*nf**2 - 2024.D0/9.D0*z3*bt0inv3*Ca**2*
     &    Cf**2*nf + 848.D0/3.D0*z3*bt0inv3*Ca**3*Cf*nf - 1408.D0/9.D0*
     &    z3*bt0inv3*Ca**4*Cf - 56*z3*bt0inv2*Ca*Cf**2*nf - 280.D0/3.D0
     &    *z3*bt0inv2*Ca**2*Cf*nf + 952.D0/3.D0*z3*bt0inv2*Ca**3*Cf - 
     &    16*z2*bt0inv4*Ca*Cf**3*nf**2 - 160.D0/3.D0*z2*bt0inv4*Ca**2*
     &    Cf**2*nf**2 - 400.D0/9.D0*z2*bt0inv4*Ca**3*Cf*nf**2 + 544.D0/
     &    3.D0*z2*bt0inv4*Ca**3*Cf**2*nf + 2720.D0/9.D0*z2*bt0inv4*
     &    Ca**4*Cf*nf )
      gnnnll = gnnnll + as**2*wminv2*zlwm * (  - 4624.D0/9.D0*z2*
     &    bt0inv4*Ca**5*Cf + 160.D0/9.D0*z2*bt0inv3*Ca*Cf**2*nf**2 + 
     &    800.D0/27.D0*z2*bt0inv3*Ca**2*Cf*nf**2 - 1072.D0/9.D0*z2*
     &    bt0inv3*Ca**2*Cf**2*nf - 8080.D0/27.D0*z2*bt0inv3*Ca**3*Cf*nf
     &     + 18224.D0/27.D0*z2*bt0inv3*Ca**4*Cf + 32.D0/3.D0*z2*bt0inv2
     &    *Cf**2*nf**2 + 160.D0/9.D0*z2*bt0inv2*Ca*Cf*nf**2 - 176.D0/3.D
     &    0*z2*bt0inv2*Ca*Cf**2*nf - 1424.D0/9.D0*z2*bt0inv2*Ca**2*Cf*
     &    nf + 2992.D0/9.D0*z2*bt0inv2*Ca**3*Cf + 16*z2*bt0inv1*Cf**2*
     &    nf + 80.D0/3.D0*z2*bt0inv1*Ca*Cf*nf - 272.D0/3.D0*z2*bt0inv1*
     &    Ca**2*Cf + 176.D0/5.D0*z2**2*bt0inv3*Ca**2*Cf**2*nf + 176.D0/
     &    3.D0*z2**2*bt0inv3*Ca**3*Cf*nf - 2992.D0/15.D0*z2**2*bt0inv3*
     &    Ca**4*Cf )
      gnnnll = gnnnll + as**2*wminv2*zlwm*zlwm2 * ( 16.D0/3.D0*bt0inv5*
     &    Cf**4*nf**3 + 80.D0/3.D0*bt0inv5*Ca*Cf**3*nf**3 + 400.D0/9.D0
     &    *bt0inv5*Ca**2*Cf**2*nf**3 - 272.D0/3.D0*bt0inv5*Ca**2*Cf**3*
     &    nf**2 + 2000.D0/81.D0*bt0inv5*Ca**3*Cf*nf**3 - 2720.D0/9.D0*
     &    bt0inv5*Ca**3*Cf**2*nf**2 - 6800.D0/27.D0*bt0inv5*Ca**4*Cf*
     &    nf**2 + 4624.D0/9.D0*bt0inv5*Ca**4*Cf**2*nf + 23120.D0/27.D0*
     &    bt0inv5*Ca**5*Cf*nf - 78608.D0/81.D0*bt0inv5*Ca**6*Cf )
      gnnnll = gnnnll + as**2*w*wminv2 * (  - 88.D0/3.D0*CA**(-2)*
     &    bt0inv3*Cf*nf**2 + 64*CA**(-2)*z3*bt0inv3*Cf*nf**2 - 4*
     &    bt0inv4*Cf**3*nf**3 + 4*bt0inv4*Cf**4*nf**2 - 422.D0/27.D0*
     &    bt0inv4*Ca*Cf**2*nf**3 + 62.D0/3.D0*bt0inv4*Ca*Cf**3*nf**2 - 
     &    1210.D0/81.D0*bt0inv4*Ca**2*Cf*nf**3 + 2452.D0/27.D0*bt0inv4*
     &    Ca**2*Cf**2*nf**2 - 68.D0/3.D0*bt0inv4*Ca**2*Cf**3*nf + 3388.D
     &    0/27.D0*bt0inv4*Ca**3*Cf*nf**2 - 5540.D0/27.D0*bt0inv4*Ca**3*
     &    Cf**2*nf - 12520.D0/27.D0*bt0inv4*Ca**4*Cf*nf + 57766.D0/81.D0
     &    *bt0inv4*Ca**5*Cf + 88.D0/9.D0*bt0inv3*Cf*nf**2 + 164.D0/243.D
     &    0*bt0inv3*Cf**2*nf**3 - 314.D0/27.D0*bt0inv3*Cf**3*nf**2 + 46
     &    *bt0inv3*Cf**4*nf + 128.D0/9.D0*bt0inv3*Ca*Cf*nf - 134.D0/243.
     &    D0*bt0inv3*Ca*Cf*nf**3 - 13798.D0/243.D0*bt0inv3*Ca*Cf**2*
     &    nf**2 - 4204.D0/27.D0*bt0inv3*Ca*Cf**3*nf - 80.D0/3.D0*
     &    bt0inv3*Ca**2*Cf - 25.D0/27.D0*bt0inv3*Ca**2*Cf*nf**2 + 97253.
     &    D0/243.D0*bt0inv3*Ca**2*Cf**2*nf - 2689.D0/81.D0*bt0inv3*
     &    Ca**3*Cf*nf )
      gnnnll = gnnnll + as**2*w*wminv2 * (  - 74437.D0/243.D0*bt0inv3*
     &    Ca**4*Cf - 224.D0/27.D0*bt0inv2*Cf**2*nf**2 - 1120.D0/81.D0*
     &    bt0inv2*Ca*Cf*nf**2 + 1616.D0/27.D0*bt0inv2*Ca*Cf**2*nf + 
     &    11888.D0/81.D0*bt0inv2*Ca**2*Cf*nf - 27472.D0/81.D0*bt0inv2*
     &    Ca**3*Cf + 1856.D0/729.D0*bt0inv1*Cf*nf**2 - 1711.D0/27.D0*
     &    bt0inv1*Cf**2*nf - 62626.D0/729.D0*bt0inv1*Ca*Cf*nf + 297029.D
     &    0/729.D0*bt0inv1*Ca**2*Cf + 192*z5*bt0inv1*Ca**2*Cf - 128.D0/
     &    9.D0*z3*Cf*nf + 704.D0/9.D0*z3*Ca*Cf - 64.D0/3.D0*z3*bt0inv3*
     &    Cf*nf**2 - 64.D0/9.D0*z3*bt0inv3*Cf**3*nf**2 - 416.D0/3.D0*z3
     &    *bt0inv3*Ca*Cf*nf + 368.D0/9.D0*z3*bt0inv3*Ca*Cf**2*nf**2 + 
     &    352.D0/9.D0*z3*bt0inv3*Ca*Cf**3*nf + 704*z3*bt0inv3*Ca**2*Cf
     &     - 416.D0/9.D0*z3*bt0inv3*Ca**2*Cf*nf**2 - 2024.D0/9.D0*z3*
     &    bt0inv3*Ca**2*Cf**2*nf + 848.D0/3.D0*z3*bt0inv3*Ca**3*Cf*nf
     &     - 1408.D0/9.D0*z3*bt0inv3*Ca**4*Cf - 56*z3*bt0inv2*Ca*Cf**2*
     &    nf - 280.D0/3.D0*z3*bt0inv2*Ca**2*Cf*nf + 952.D0/3.D0*z3*
     &    bt0inv2*Ca**3*Cf )
      gnnnll = gnnnll + as**2*w*wminv2 * (  - 160.D0/27.D0*z3*bt0inv1*
     &    Cf*nf**2 + 304.D0/9.D0*z3*bt0inv1*Cf**2*nf + 1240.D0/9.D0*z3*
     &    bt0inv1*Ca*Cf*nf - 20072.D0/27.D0*z3*bt0inv1*Ca**2*Cf - 160.D0
     &    /9.D0*z2*Cf*nf + 1072.D0/9.D0*z2*Ca*Cf - 16*z2*bt0inv4*Ca*
     &    Cf**3*nf**2 - 160.D0/3.D0*z2*bt0inv4*Ca**2*Cf**2*nf**2 - 400.D
     &    0/9.D0*z2*bt0inv4*Ca**3*Cf*nf**2 + 544.D0/3.D0*z2*bt0inv4*
     &    Ca**3*Cf**2*nf + 2720.D0/9.D0*z2*bt0inv4*Ca**4*Cf*nf - 4624.D0
     &    /9.D0*z2*bt0inv4*Ca**5*Cf + 160.D0/9.D0*z2*bt0inv3*Ca*Cf**2*
     &    nf**2 + 800.D0/27.D0*z2*bt0inv3*Ca**2*Cf*nf**2 - 1072.D0/9.D0
     &    *z2*bt0inv3*Ca**2*Cf**2*nf - 8080.D0/27.D0*z2*bt0inv3*Ca**3*
     &    Cf*nf + 18224.D0/27.D0*z2*bt0inv3*Ca**4*Cf + 32.D0/3.D0*z2*
     &    bt0inv2*Cf**2*nf**2 + 160.D0/9.D0*z2*bt0inv2*Ca*Cf*nf**2 - 
     &    176.D0/3.D0*z2*bt0inv2*Ca*Cf**2*nf - 1424.D0/9.D0*z2*bt0inv2*
     &    Ca**2*Cf*nf + 2992.D0/9.D0*z2*bt0inv2*Ca**3*Cf - 320.D0/27.D0
     &    *z2*bt0inv1*Cf*nf**2 + 16*z2*bt0inv1*Cf**2*nf + 14696.D0/81.D0
     &    *z2*bt0inv1*Ca*Cf*nf )
      gnnnll = gnnnll + as**2*w*wminv2 * (  - 49112.D0/81.D0*z2*bt0inv1
     &    *Ca**2*Cf + 176.D0/3.D0*z2*z3*bt0inv1*Ca**2*Cf - 32*z2**2*Ca*
     &    Cf + 176.D0/5.D0*z2**2*bt0inv3*Ca**2*Cf**2*nf + 176.D0/3.D0*
     &    z2**2*bt0inv3*Ca**3*Cf*nf - 2992.D0/15.D0*z2**2*bt0inv3*Ca**4
     &    *Cf + 32.D0/5.D0*z2**2*bt0inv1*Cf**2*nf - 368.D0/15.D0*z2**2*
     &    bt0inv1*Ca*Cf*nf + 1496.D0/15.D0*z2**2*bt0inv1*Ca**2*Cf )
      gnnnll = gnnnll + as**2*w*wminv2*zlwm * ( 176.D0/3.D0*CA**(-2)*
     &    bt0inv3*Cf*nf**2 - 128*CA**(-2)*z3*bt0inv3*Cf*nf**2 - 88.D0/9.
     &    D0*bt0inv4*Cf**3*nf**3 - 8*bt0inv4*Cf**4*nf**2 - 28*bt0inv4*
     &    Ca*Cf**2*nf**3 + 700.D0/9.D0*bt0inv4*Ca*Cf**3*nf**2 - 1580.D0/
     &    81.D0*bt0inv4*Ca**2*Cf*nf**3 + 3752.D0/9.D0*bt0inv4*Ca**2*
     &    Cf**2*nf**2 + 136.D0/3.D0*bt0inv4*Ca**2*Cf**3*nf + 11224.D0/
     &    27.D0*bt0inv4*Ca**3*Cf*nf**2 - 8456.D0/9.D0*bt0inv4*Ca**3*
     &    Cf**2*nf - 5680.D0/3.D0*bt0inv4*Ca**4*Cf*nf + 194276.D0/81.D0
     &    *bt0inv4*Ca**5*Cf - 176.D0/9.D0*bt0inv3*Cf*nf**2 - 616.D0/243.
     &    D0*bt0inv3*Cf**2*nf**3 - 1352.D0/27.D0*bt0inv3*Cf**3*nf**2 - 
     &    92*bt0inv3*Cf**4*nf - 256.D0/9.D0*bt0inv3*Ca*Cf*nf - 212.D0/
     &    243.D0*bt0inv3*Ca*Cf*nf**3 - 17152.D0/243.D0*bt0inv3*Ca*Cf**2
     &    *nf**2 + 8408.D0/27.D0*bt0inv3*Ca*Cf**3*nf + 160.D0/3.D0*
     &    bt0inv3*Ca**2*Cf - 7666.D0/81.D0*bt0inv3*Ca**2*Cf*nf**2 - 
     &    14146.D0/243.D0*bt0inv3*Ca**2*Cf**2*nf + 77902.D0/81.D0*
     &    bt0inv3*Ca**3*Cf*nf )
      gnnnll = gnnnll + as**2*w*wminv2*zlwm * (  - 300946.D0/243.D0*
     &    bt0inv3*Ca**4*Cf + 128.D0/3.D0*z3*bt0inv3*Cf*nf**2 + 704.D0/9.
     &    D0*z3*bt0inv3*Cf**3*nf**2 + 832.D0/3.D0*z3*bt0inv3*Ca*Cf*nf
     &     - 448.D0/9.D0*z3*bt0inv3*Ca*Cf**2*nf**2 - 704.D0/9.D0*z3*
     &    bt0inv3*Ca*Cf**3*nf - 1408*z3*bt0inv3*Ca**2*Cf - 32*z3*
     &    bt0inv3*Ca**2*Cf*nf**2 + 1312.D0/9.D0*z3*bt0inv3*Ca**2*Cf**2*
     &    nf - 400.D0/9.D0*z3*bt0inv3*Ca**3*Cf*nf - 176.D0/9.D0*z3*
     &    bt0inv3*Ca**4*Cf )
      gnnnll = gnnnll + as**2*w**2*wminv2 * ( 44*CA**(-2)*bt0inv3*Cf*
     &    nf**2 - 96*CA**(-2)*z3*bt0inv3*Cf*nf**2 - 26.D0/9.D0*bt0inv4*
     &    Cf**3*nf**3 - 6*bt0inv4*Cf**4*nf**2 - 167.D0/27.D0*bt0inv4*Ca
     &    *Cf**2*nf**3 + 257.D0/9.D0*bt0inv4*Ca*Cf**3*nf**2 - 185.D0/81.
     &    D0*bt0inv4*Ca**2*Cf*nf**3 + 4402.D0/27.D0*bt0inv4*Ca**2*Cf**2
     &    *nf**2 + 34*bt0inv4*Ca**2*Cf**3*nf + 1306.D0/9.D0*bt0inv4*
     &    Ca**3*Cf*nf**2 - 9914.D0/27.D0*bt0inv4*Ca**3*Cf**2*nf - 19300.
     &    D0/27.D0*bt0inv4*Ca**4*Cf*nf + 68255.D0/81.D0*bt0inv4*Ca**5*
     &    Cf - 44.D0/3.D0*bt0inv3*Cf*nf**2 - 178.D0/81.D0*bt0inv3*Cf**2
     &    *nf**3 - 503.D0/9.D0*bt0inv3*Cf**3*nf**2 - 69*bt0inv3*Cf**4*
     &    nf - 64.D0/3.D0*bt0inv3*Ca*Cf*nf - 31.D0/27.D0*bt0inv3*Ca*Cf*
     &    nf**3 - 8017.D0/81.D0*bt0inv3*Ca*Cf**2*nf**2 + 2102.D0/9.D0*
     &    bt0inv3*Ca*Cf**3*nf + 40*bt0inv3*Ca**2*Cf - 15407.D0/162.D0*
     &    bt0inv3*Ca**2*Cf*nf**2 + 22987.D0/162.D0*bt0inv3*Ca**2*Cf**2*
     &    nf + 153115.D0/162.D0*bt0inv3*Ca**3*Cf*nf - 225443.D0/162.D0*
     &    bt0inv3*Ca**4*Cf )
      gnnnll = gnnnll + as**2*w**2*wminv2 * ( 10351*bt0inv2 - 51719.D0/
     &    20.D0*bt0inv2*nf + 488943.D0/5000.D0*bt0inv2*nf**2 + 409043.D0
     &    /250000.D0*bt0inv2*nf**3 + 112.D0/27.D0*bt0inv2*Cf**2*nf**2
     &     + 560.D0/81.D0*bt0inv2*Ca*Cf*nf**2 - 808.D0/27.D0*bt0inv2*Ca
     &    *Cf**2*nf - 5944.D0/81.D0*bt0inv2*Ca**2*Cf*nf + 13736.D0/81.D0
     &    *bt0inv2*Ca**3*Cf - 928.D0/729.D0*bt0inv1*Cf*nf**2 + 1711.D0/
     &    54.D0*bt0inv1*Cf**2*nf + 31313.D0/729.D0*bt0inv1*Ca*Cf*nf - 
     &    297029.D0/1458.D0*bt0inv1*Ca**2*Cf - 96*z5*bt0inv1*Ca**2*Cf
     &     + 64.D0/9.D0*z3*Cf*nf - 352.D0/9.D0*z3*Ca*Cf + 32*z3*bt0inv3
     &    *Cf*nf**2 + 224.D0/3.D0*z3*bt0inv3*Cf**3*nf**2 + 208*z3*
     &    bt0inv3*Ca*Cf*nf - 88.D0/3.D0*z3*bt0inv3*Ca*Cf**2*nf**2 - 176.
     &    D0/3.D0*z3*bt0inv3*Ca*Cf**3*nf - 1056*z3*bt0inv3*Ca**2*Cf - 
     &    496.D0/9.D0*z3*bt0inv3*Ca**2*Cf*nf**2 + 100.D0/3.D0*z3*
     &    bt0inv3*Ca**2*Cf**2*nf + 872.D0/9.D0*z3*bt0inv3*Ca**3*Cf*nf
     &     - 880.D0/9.D0*z3*bt0inv3*Ca**4*Cf + 28*z3*bt0inv2*Ca*Cf**2*
     &    nf )
      gnnnll = gnnnll + as**2*w**2*wminv2 * ( 140.D0/3.D0*z3*bt0inv2*
     &    Ca**2*Cf*nf - 476.D0/3.D0*z3*bt0inv2*Ca**3*Cf + 80.D0/27.D0*
     &    z3*bt0inv1*Cf*nf**2 - 152.D0/9.D0*z3*bt0inv1*Cf**2*nf - 620.D0
     &    /9.D0*z3*bt0inv1*Ca*Cf*nf + 10036.D0/27.D0*z3*bt0inv1*Ca**2*
     &    Cf + 80.D0/9.D0*z2*Cf*nf - 536.D0/9.D0*z2*Ca*Cf + 8*z2*
     &    bt0inv4*Ca*Cf**3*nf**2 + 80.D0/3.D0*z2*bt0inv4*Ca**2*Cf**2*
     &    nf**2 + 200.D0/9.D0*z2*bt0inv4*Ca**3*Cf*nf**2 - 272.D0/3.D0*
     &    z2*bt0inv4*Ca**3*Cf**2*nf - 1360.D0/9.D0*z2*bt0inv4*Ca**4*Cf*
     &    nf + 2312.D0/9.D0*z2*bt0inv4*Ca**5*Cf + 80.D0/9.D0*z2*bt0inv3
     &    *Ca*Cf**2*nf**2 + 400.D0/27.D0*z2*bt0inv3*Ca**2*Cf*nf**2 - 
     &    536.D0/9.D0*z2*bt0inv3*Ca**2*Cf**2*nf - 4040.D0/27.D0*z2*
     &    bt0inv3*Ca**3*Cf*nf + 9112.D0/27.D0*z2*bt0inv3*Ca**4*Cf - 16.D
     &    0/3.D0*z2*bt0inv2*Cf**2*nf**2 - 80.D0/9.D0*z2*bt0inv2*Ca*Cf*
     &    nf**2 + 88.D0/3.D0*z2*bt0inv2*Ca*Cf**2*nf + 712.D0/9.D0*z2*
     &    bt0inv2*Ca**2*Cf*nf - 1496.D0/9.D0*z2*bt0inv2*Ca**3*Cf + 160.D
     &    0/27.D0*z2*bt0inv1*Cf*nf**2 )
      gnnnll = gnnnll + as**2*w**2*wminv2 * (  - 8*z2*bt0inv1*Cf**2*nf
     &     - 7348.D0/81.D0*z2*bt0inv1*Ca*Cf*nf + 24556.D0/81.D0*z2*
     &    bt0inv1*Ca**2*Cf - 88.D0/3.D0*z2*z3*bt0inv1*Ca**2*Cf + 16*
     &    z2**2*Ca*Cf + 88.D0/5.D0*z2**2*bt0inv3*Ca**2*Cf**2*nf + 88.D0/
     &    3.D0*z2**2*bt0inv3*Ca**3*Cf*nf - 1496.D0/15.D0*z2**2*bt0inv3*
     &    Ca**4*Cf - 16.D0/5.D0*z2**2*bt0inv1*Cf**2*nf + 184.D0/15.D0*
     &    z2**2*bt0inv1*Ca*Cf*nf - 748.D0/15.D0*z2**2*bt0inv1*Ca**2*Cf
     &     )
      gnnnll = gnnnll + as**2*w**2*wminv2*zlwm * (  - 88.D0/3.D0*
     &    CA**(-2)*bt0inv3*Cf*nf**2 + 64*CA**(-2)*z3*bt0inv3*Cf*nf**2
     &     - 16*bt0inv5*Cf**4*nf**3 - 80*bt0inv5*Ca*Cf**3*nf**3 - 400.D0
     &    /3.D0*bt0inv5*Ca**2*Cf**2*nf**3 + 272*bt0inv5*Ca**2*Cf**3*
     &    nf**2 - 2000.D0/27.D0*bt0inv5*Ca**3*Cf*nf**3 + 2720.D0/3.D0*
     &    bt0inv5*Ca**3*Cf**2*nf**2 + 6800.D0/9.D0*bt0inv5*Ca**4*Cf*
     &    nf**2 - 4624.D0/3.D0*bt0inv5*Ca**4*Cf**2*nf - 23120.D0/9.D0*
     &    bt0inv5*Ca**5*Cf*nf + 78608.D0/27.D0*bt0inv5*Ca**6*Cf + 88.D0/
     &    9.D0*bt0inv4*Cf**3*nf**3 + 8*bt0inv4*Cf**4*nf**2 + 28*bt0inv4
     &    *Ca*Cf**2*nf**3 - 700.D0/9.D0*bt0inv4*Ca*Cf**3*nf**2 + 1580.D0
     &    /81.D0*bt0inv4*Ca**2*Cf*nf**3 - 3752.D0/9.D0*bt0inv4*Ca**2*
     &    Cf**2*nf**2 - 136.D0/3.D0*bt0inv4*Ca**2*Cf**3*nf - 11224.D0/
     &    27.D0*bt0inv4*Ca**3*Cf*nf**2 + 8456.D0/9.D0*bt0inv4*Ca**3*
     &    Cf**2*nf + 5680.D0/3.D0*bt0inv4*Ca**4*Cf*nf - 194276.D0/81.D0
     &    *bt0inv4*Ca**5*Cf + 88.D0/9.D0*bt0inv3*Cf*nf**2 + 308.D0/243.D
     &    0*bt0inv3*Cf**2*nf**3 )
      gnnnll = gnnnll + as**2*w**2*wminv2*zlwm * ( 676.D0/27.D0*bt0inv3
     &    *Cf**3*nf**2 + 46*bt0inv3*Cf**4*nf + 128.D0/9.D0*bt0inv3*Ca*
     &    Cf*nf + 106.D0/243.D0*bt0inv3*Ca*Cf*nf**3 + 8576.D0/243.D0*
     &    bt0inv3*Ca*Cf**2*nf**2 - 4204.D0/27.D0*bt0inv3*Ca*Cf**3*nf - 
     &    80.D0/3.D0*bt0inv3*Ca**2*Cf + 3833.D0/81.D0*bt0inv3*Ca**2*Cf*
     &    nf**2 + 7073.D0/243.D0*bt0inv3*Ca**2*Cf**2*nf - 38951.D0/81.D0
     &    *bt0inv3*Ca**3*Cf*nf + 150473.D0/243.D0*bt0inv3*Ca**4*Cf - 64.
     &    D0/3.D0*z3*bt0inv3*Cf*nf**2 - 352.D0/9.D0*z3*bt0inv3*Cf**3*
     &    nf**2 - 416.D0/3.D0*z3*bt0inv3*Ca*Cf*nf + 224.D0/9.D0*z3*
     &    bt0inv3*Ca*Cf**2*nf**2 + 352.D0/9.D0*z3*bt0inv3*Ca*Cf**3*nf
     &     + 704*z3*bt0inv3*Ca**2*Cf + 16*z3*bt0inv3*Ca**2*Cf*nf**2 - 
     &    656.D0/9.D0*z3*bt0inv3*Ca**2*Cf**2*nf + 200.D0/9.D0*z3*
     &    bt0inv3*Ca**3*Cf*nf + 88.D0/9.D0*z3*bt0inv3*Ca**4*Cf )
      gnnnll = gnnnll + as**2*w**3*wminv2 * (  - 176.D0/9.D0*CA**(-2)*
     &    bt0inv3*Cf*nf**2 + 128.D0/3.D0*CA**(-2)*z3*bt0inv3*Cf*nf**2
     &     - 32.D0/3.D0*bt0inv5*Cf**4*nf**3 - 160.D0/3.D0*bt0inv5*Ca*
     &    Cf**3*nf**3 - 800.D0/9.D0*bt0inv5*Ca**2*Cf**2*nf**3 + 544.D0/
     &    3.D0*bt0inv5*Ca**2*Cf**3*nf**2 - 4000.D0/81.D0*bt0inv5*Ca**3*
     &    Cf*nf**3 + 5440.D0/9.D0*bt0inv5*Ca**3*Cf**2*nf**2 + 13600.D0/
     &    27.D0*bt0inv5*Ca**4*Cf*nf**2 - 9248.D0/9.D0*bt0inv5*Ca**4*
     &    Cf**2*nf - 46240.D0/27.D0*bt0inv5*Ca**5*Cf*nf + 157216.D0/81.D
     &    0*bt0inv5*Ca**6*Cf + 112.D0/9.D0*bt0inv4*Cf**3*nf**3 + 16.D0/
     &    3.D0*bt0inv4*Cf**4*nf**2 + 3112.D0/81.D0*bt0inv4*Ca*Cf**2*
     &    nf**3 - 824.D0/9.D0*bt0inv4*Ca*Cf**3*nf**2 + 7160.D0/243.D0*
     &    bt0inv4*Ca**2*Cf*nf**3 - 38672.D0/81.D0*bt0inv4*Ca**2*Cf**2*
     &    nf**2 - 272.D0/9.D0*bt0inv4*Ca**2*Cf**3*nf - 40448.D0/81.D0*
     &    bt0inv4*Ca**3*Cf*nf**2 + 87184.D0/81.D0*bt0inv4*Ca**3*Cf**2*
     &    nf + 178400.D0/81.D0*bt0inv4*Ca**4*Cf*nf - 698360.D0/243.D0*
     &    bt0inv4*Ca**5*Cf )
      gnnnll = gnnnll + as**2*w**3*wminv2 * ( 176.D0/27.D0*bt0inv3*Cf*
     &    nf**2 - 416.D0/729.D0*bt0inv3*Cf**2*nf**3 + 3212.D0/81.D0*
     &    bt0inv3*Cf**3*nf**2 + 92.D0/3.D0*bt0inv3*Cf**4*nf + 256.D0/27.
     &    D0*bt0inv3*Ca*Cf*nf - 296.D0/243.D0*bt0inv3*Ca*Cf*nf**3 + 
     &    83044.D0/729.D0*bt0inv3*Ca*Cf**2*nf**2 - 7604.D0/81.D0*
     &    bt0inv3*Ca*Cf**3*nf - 160.D0/9.D0*bt0inv3*Ca**2*Cf + 28444.D0/
     &    243.D0*bt0inv3*Ca**2*Cf*nf**2 - 248624.D0/729.D0*bt0inv3*
     &    Ca**2*Cf**2*nf - 232676.D0/243.D0*bt0inv3*Ca**3*Cf*nf + 
     &    125956.D0/81.D0*bt0inv3*Ca**4*Cf - 20702.D0/3.D0*bt0inv2 + 
     &    51719.D0/30.D0*bt0inv2*nf - 162981.D0/2500.D0*bt0inv2*nf**2
     &     - 409043.D0/375000.D0*bt0inv2*nf**3 - 128.D0/9.D0*z3*bt0inv3
     &    *Cf*nf**2 - 1280.D0/27.D0*z3*bt0inv3*Cf**3*nf**2 - 832.D0/9.D0
     &    *z3*bt0inv3*Ca*Cf*nf + 160.D0/27.D0*z3*bt0inv3*Ca*Cf**2*nf**2
     &     + 704.D0/27.D0*z3*bt0inv3*Ca*Cf**3*nf + 1408.D0/3.D0*z3*
     &    bt0inv3*Ca**2*Cf + 1408.D0/27.D0*z3*bt0inv3*Ca**2*Cf*nf**2 + 
     &    1424.D0/27.D0*z3*bt0inv3*Ca**2*Cf**2*nf )
      gnnnll = gnnnll + as**2*w**3*wminv2 * (  - 4288.D0/27.D0*z3*
     &    bt0inv3*Ca**3*Cf*nf + 352.D0/3.D0*z3*bt0inv3*Ca**4*Cf + 32.D0/
     &    3.D0*z2*bt0inv4*Ca*Cf**3*nf**2 + 320.D0/9.D0*z2*bt0inv4*Ca**2
     &    *Cf**2*nf**2 + 800.D0/27.D0*z2*bt0inv4*Ca**3*Cf*nf**2 - 1088.D
     &    0/9.D0*z2*bt0inv4*Ca**3*Cf**2*nf - 5440.D0/27.D0*z2*bt0inv4*
     &    Ca**4*Cf*nf + 9248.D0/27.D0*z2*bt0inv4*Ca**5*Cf - 136.D0/9.D0
     &    *z2*bt0inv3*Ca*Cf**2*nf**2 - 8.D0/3.D0*z2*bt0inv3*Ca*Cf**3*nf
     &     - 1916.D0/81.D0*z2*bt0inv3*Ca**2*Cf*nf**2 + 988.D0/9.D0*z2*
     &    bt0inv3*Ca**2*Cf**2*nf + 21820.D0/81.D0*z2*bt0inv3*Ca**3*Cf*
     &    nf - 47876.D0/81.D0*z2*bt0inv3*Ca**4*Cf - 352.D0/15.D0*z2**2*
     &    bt0inv3*Ca**2*Cf**2*nf - 352.D0/9.D0*z2**2*bt0inv3*Ca**3*Cf*
     &    nf + 5984.D0/45.D0*z2**2*bt0inv3*Ca**4*Cf )
      gnnnll = gnnnll + lqr*zlwm * ( 4*bt0inv1*Cf )
      gnnnll = gnnnll + lqr*as*wminv*zlwm * (  - 8*bt0inv2*Cf**2*nf - 
     &    40.D0/3.D0*bt0inv2*Ca*Cf*nf + 136.D0/3.D0*bt0inv2*Ca**2*Cf )
      gnnnll = gnnnll + lqr*as*w*wminv * (  - 8*bt0inv2*Cf**2*nf - 40.D0
     &    /3.D0*bt0inv2*Ca*Cf*nf + 136.D0/3.D0*bt0inv2*Ca**2*Cf + 40.D0/
     &    9.D0*bt0inv1*Cf*nf - 268.D0/9.D0*bt0inv1*Ca*Cf + 8*z2*bt0inv1
     &    *Ca*Cf )
      gnnnll = gnnnll + lqr*as**2*wminv2*zlwm2 * (  - 8*bt0inv3*Cf**3*
     &    nf**2 - 80.D0/3.D0*bt0inv3*Ca*Cf**2*nf**2 - 200.D0/9.D0*
     &    bt0inv3*Ca**2*Cf*nf**2 + 272.D0/3.D0*bt0inv3*Ca**2*Cf**2*nf
     &     + 1360.D0/9.D0*bt0inv3*Ca**3*Cf*nf - 2312.D0/9.D0*bt0inv3*
     &    Ca**4*Cf )
      gnnnll = gnnnll + lqr*as**2*wminv2*zlwm * ( 80.D0/9.D0*bt0inv2*
     &    Cf**2*nf**2 + 400.D0/27.D0*bt0inv2*Ca*Cf*nf**2 - 536.D0/9.D0*
     &    bt0inv2*Ca*Cf**2*nf - 4040.D0/27.D0*bt0inv2*Ca**2*Cf*nf + 
     &    9112.D0/27.D0*bt0inv2*Ca**3*Cf + 16*z2*bt0inv2*Ca*Cf**2*nf + 
     &    80.D0/3.D0*z2*bt0inv2*Ca**2*Cf*nf - 272.D0/3.D0*z2*bt0inv2*
     &    Ca**3*Cf )
      gnnnll = gnnnll + lqr*as**2*w*wminv2 * ( 224.D0/27.D0*Cf*nf - 
     &    1616.D0/27.D0*Ca*Cf + 80.D0/9.D0*bt0inv2*Cf**2*nf**2 + 400.D0/
     &    27.D0*bt0inv2*Ca*Cf*nf**2 - 536.D0/9.D0*bt0inv2*Ca*Cf**2*nf
     &     - 4040.D0/27.D0*bt0inv2*Ca**2*Cf*nf + 9112.D0/27.D0*bt0inv2*
     &    Ca**3*Cf + 16.D0/27.D0*bt0inv1*Cf*nf**2 + 110.D0/3.D0*bt0inv1
     &    *Cf**2*nf + 836.D0/27.D0*bt0inv1*Ca*Cf*nf - 490.D0/3.D0*
     &    bt0inv1*Ca**2*Cf + 56*z3*Ca*Cf - 32*z3*bt0inv1*Cf**2*nf + 112.
     &    D0/3.D0*z3*bt0inv1*Ca*Cf*nf - 88.D0/3.D0*z3*bt0inv1*Ca**2*Cf
     &     + 16*z2*bt0inv2*Ca*Cf**2*nf + 80.D0/3.D0*z2*bt0inv2*Ca**2*Cf
     &    *nf - 272.D0/3.D0*z2*bt0inv2*Ca**3*Cf - 160.D0/9.D0*z2*
     &    bt0inv1*Ca*Cf*nf + 1072.D0/9.D0*z2*bt0inv1*Ca**2*Cf - 176.D0/
     &    5.D0*z2**2*bt0inv1*Ca**2*Cf )
      gnnnll = gnnnll + lqr*as**2*w**2*wminv2 * (  - 112.D0/27.D0*Cf*nf
     &     + 808.D0/27.D0*Ca*Cf + 8*bt0inv3*Cf**3*nf**2 + 80.D0/3.D0*
     &    bt0inv3*Ca*Cf**2*nf**2 + 200.D0/9.D0*bt0inv3*Ca**2*Cf*nf**2
     &     - 272.D0/3.D0*bt0inv3*Ca**2*Cf**2*nf - 1360.D0/9.D0*bt0inv3*
     &    Ca**3*Cf*nf + 2312.D0/9.D0*bt0inv3*Ca**4*Cf - 62.D0/9.D0*
     &    bt0inv2*Cf**2*nf**2 - 2*bt0inv2*Cf**3*nf - 31.D0/3.D0*bt0inv2
     &    *Ca*Cf*nf**2 + 473.D0/9.D0*bt0inv2*Ca*Cf**2*nf + 1145.D0/9.D0
     &    *bt0inv2*Ca**2*Cf*nf - 2471.D0/9.D0*bt0inv2*Ca**3*Cf - 8.D0/
     &    27.D0*bt0inv1*Cf*nf**2 - 55.D0/3.D0*bt0inv1*Cf**2*nf - 418.D0/
     &    27.D0*bt0inv1*Ca*Cf*nf + 245.D0/3.D0*bt0inv1*Ca**2*Cf - 28*z3
     &    *Ca*Cf + 16*z3*bt0inv1*Cf**2*nf - 56.D0/3.D0*z3*bt0inv1*Ca*Cf
     &    *nf + 44.D0/3.D0*z3*bt0inv1*Ca**2*Cf - 8*z2*bt0inv2*Ca*Cf**2*
     &    nf - 40.D0/3.D0*z2*bt0inv2*Ca**2*Cf*nf + 136.D0/3.D0*z2*
     &    bt0inv2*Ca**3*Cf + 80.D0/9.D0*z2*bt0inv1*Ca*Cf*nf - 536.D0/9.D
     &    0*z2*bt0inv1*Ca**2*Cf + 88.D0/5.D0*z2**2*bt0inv1*Ca**2*Cf )
      gnnnll = gnnnll + lqr**2*as*w*wminv * ( 2*Cf )
      gnnnll = gnnnll + lqr**2*as**2*wminv2*zlwm * ( 4*bt0inv1*Cf**2*nf
     &     + 20.D0/3.D0*bt0inv1*Ca*Cf*nf - 68.D0/3.D0*bt0inv1*Ca**2*Cf
     &     )
      gnnnll = gnnnll + lqr**2*as**2*w*wminv2 * (  - 40.D0/9.D0*Cf*nf
     &     + 268.D0/9.D0*Ca*Cf - 8*z2*Ca*Cf )
      gnnnll = gnnnll + lqr**2*as**2*w**2*wminv2 * ( 20.D0/9.D0*Cf*nf
     &     - 134.D0/9.D0*Ca*Cf + 4*z2*Ca*Cf )
      gnnnll = gnnnll + lqr**3*as**2*w*wminv2 * ( 8.D0/9.D0*Cf*nf - 44.D
     &    0/9.D0*Ca*Cf )
      gnnnll = gnnnll + lqr**3*as**2*w**2*wminv2 * (  - 4.D0/9.D0*Cf*nf
     &     + 22.D0/9.D0*Ca*Cf )
      gnnnll = gnnnll + lfr*w * ( 4*bt0inv1*Cf )
      gnnnll = gnnnll + lfr*as*w*wminv * (  - 40.D0/9.D0*bt0inv1*Cf*nf
     &     + 268.D0/9.D0*bt0inv1*Ca*Cf - 8*z2*bt0inv1*Ca*Cf )
      gnnnll = gnnnll + lfr*as*w**2*wminv * ( 40.D0/9.D0*bt0inv1*Cf*nf
     &     - 268.D0/9.D0*bt0inv1*Ca*Cf + 8*z2*bt0inv1*Ca*Cf )
      gnnnll = gnnnll + lfr*as**2*w*wminv2 * (  - 16.D0/27.D0*bt0inv1*
     &    Cf*nf**2 - 110.D0/3.D0*bt0inv1*Cf**2*nf - 836.D0/27.D0*
     &    bt0inv1*Ca*Cf*nf + 490.D0/3.D0*bt0inv1*Ca**2*Cf + 32*z3*
     &    bt0inv1*Cf**2*nf - 112.D0/3.D0*z3*bt0inv1*Ca*Cf*nf + 88.D0/3.D
     &    0*z3*bt0inv1*Ca**2*Cf + 160.D0/9.D0*z2*bt0inv1*Ca*Cf*nf - 
     &    1072.D0/9.D0*z2*bt0inv1*Ca**2*Cf + 176.D0/5.D0*z2**2*bt0inv1*
     &    Ca**2*Cf )
      gnnnll = gnnnll + lfr*as**2*w**2*wminv2 * ( 32.D0/27.D0*bt0inv1*
     &    Cf*nf**2 + 220.D0/3.D0*bt0inv1*Cf**2*nf + 1672.D0/27.D0*
     &    bt0inv1*Ca*Cf*nf - 980.D0/3.D0*bt0inv1*Ca**2*Cf - 64*z3*
     &    bt0inv1*Cf**2*nf + 224.D0/3.D0*z3*bt0inv1*Ca*Cf*nf - 176.D0/3.
     &    D0*z3*bt0inv1*Ca**2*Cf - 320.D0/9.D0*z2*bt0inv1*Ca*Cf*nf + 
     &    2144.D0/9.D0*z2*bt0inv1*Ca**2*Cf - 352.D0/5.D0*z2**2*bt0inv1*
     &    Ca**2*Cf )
      gnnnll = gnnnll + lfr*as**2*w**3*wminv2 * (  - 16.D0/27.D0*
     &    bt0inv1*Cf*nf**2 - 110.D0/3.D0*bt0inv1*Cf**2*nf - 836.D0/27.D0
     &    *bt0inv1*Ca*Cf*nf + 490.D0/3.D0*bt0inv1*Ca**2*Cf + 32*z3*
     &    bt0inv1*Cf**2*nf - 112.D0/3.D0*z3*bt0inv1*Ca*Cf*nf + 88.D0/3.D
     &    0*z3*bt0inv1*Ca**2*Cf + 160.D0/9.D0*z2*bt0inv1*Ca*Cf*nf - 
     &    1072.D0/9.D0*z2*bt0inv1*Ca**2*Cf + 176.D0/5.D0*z2**2*bt0inv1*
     &    Ca**2*Cf )
      gnnnll = gnnnll + lfr**2*as*w*wminv * (  - 2*Cf )
      gnnnll = gnnnll + lfr**2*as*w**2*wminv * ( 2*Cf )
      gnnnll = gnnnll + lfr**2*as**2*w*wminv2 * ( 40.D0/9.D0*Cf*nf - 
     &    268.D0/9.D0*Ca*Cf + 4*bt0inv1*Cf**2*nf + 20.D0/3.D0*bt0inv1*
     &    Ca*Cf*nf - 68.D0/3.D0*bt0inv1*Ca**2*Cf + 8*z2*Ca*Cf )
      gnnnll = gnnnll + lfr**2*as**2*w**2*wminv2 * (  - 80.D0/9.D0*Cf*
     &    nf + 536.D0/9.D0*Ca*Cf - 8*bt0inv1*Cf**2*nf - 40.D0/3.D0*
     &    bt0inv1*Ca*Cf*nf + 136.D0/3.D0*bt0inv1*Ca**2*Cf - 16*z2*Ca*Cf
     &     )
      gnnnll = gnnnll + lfr**2*as**2*w**3*wminv2 * ( 40.D0/9.D0*Cf*nf
     &     - 268.D0/9.D0*Ca*Cf + 4*bt0inv1*Cf**2*nf + 20.D0/3.D0*
     &    bt0inv1*Ca*Cf*nf - 68.D0/3.D0*bt0inv1*Ca**2*Cf + 8*z2*Ca*Cf )
      gnnnll = gnnnll + lfr**3*as**2*w*wminv2 * (  - 8.D0/9.D0*Cf*nf + 
     &    44.D0/9.D0*Ca*Cf )
      gnnnll = gnnnll + lfr**3*as**2*w**2*wminv2 * ( 16.D0/9.D0*Cf*nf
     &     - 88.D0/9.D0*Ca*Cf )
      gnnnll = gnnnll + lfr**3*as**2*w**3*wminv2 * (  - 8.D0/9.D0*Cf*nf
     &     + 44.D0/9.D0*Ca*Cf )
      gnnnll = gnnnll + lnNb * ( 8*bt0inv1*Cf )
      gnnnll = gnnnll + lnNb*zlwm * (  - 8*bt0inv1*Cf )
      gnnnll = gnnnll + lnNb*winv*zlwm * ( 8*bt0inv1*Cf ) 

       CHg0tnnnll1=  + as * (  - 16*Cf + 16*z2*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + as**2 * ( 127.D0/6.D0*Cf*nf + 511.D0/
     &    4.D0*Cf**2 - 1535.D0/12.D0*Ca*Cf + 8.D0/9.D0*z3*Cf*nf - 60*z3
     &    *Cf**2 + 604.D0/9.D0*z3*Ca*Cf - 64.D0/3.D0*z2*Cf*nf - 198*z2*
     &    Cf**2 + 376.D0/3.D0*z2*Ca*Cf + 552.D0/5.D0*z2**2*Cf**2 - 92.D0
     &    /5.D0*z2**2*Ca*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + as**3 * (  - 7081.D0/243.D0*Cf*nf**2
     &     - 421.D0/3.D0*Cf**2*nf - 5599.D0/6.D0*Cf**3 + 110651.D0/243.D
     &    0*Ca*Cf*nf + 74321.D0/36.D0*Ca*Cf**2 - 1505881.D0/972.D0*
     &    Ca**2*Cf - 608.D0/9.D0*z5*Cf**2*nf + 1328*z5*Cf**3 - 8*z5*Ca*
     &    Cf*nf - 5512.D0/9.D0*z5*Ca*Cf**2 - 204*z5*Ca**2*Cf + 16.D0/81.
     &    D0*z3*Cf*nf**2 + 9448.D0/27.D0*z3*Cf**2*nf - 460*z3*Cf**3 - 
     &    24512.D0/81.D0*z3*Ca*Cf*nf - 51508.D0/27.D0*z3*Ca*Cf**2 + 
     &    139345.D0/81.D0*z3*Ca**2*Cf + 32*z3**2*Cf**3 + 592.D0/3.D0*
     &    z3**2*Ca*Cf**2 - 400.D0/3.D0*z3**2*Ca**2*Cf + 1072.D0/27.D0*
     &    z2*Cf*nf**2 + 9064.D0/27.D0*z2*Cf**2*nf + 2936.D0/3.D0*z2*
     &    Cf**3 - 44540.D0/81.D0*z2*Ca*Cf*nf - 66544.D0/27.D0*z2*Ca*
     &    Cf**2 + 130295.D0/81.D0*z2*Ca**2*Cf - 256.D0/3.D0*z2*z3*Cf**2
     &    *nf - 400*z2*z3*Cf**3 + 880.D0/9.D0*z2*z3*Ca*Cf*nf + 3680.D0/
     &    3.D0*z2*z3*Ca*Cf**2 - 7228.D0/9.D0*z2*z3*Ca**2*Cf + 448.D0/
     &    135.D0*z2**2*Cf*nf**2 - 36208.D0/135.D0*z2**2*Cf**2*nf - 5972.
     &    D0/5.D0*z2**2*Cf**3 )
      CHg0tnnnll1 = CHg0tnnnll1 + as**3 * ( 1156.D0/135.D0*z2**2*Ca*Cf*
     &    nf + 258304.D0/135.D0*z2**2*Ca*Cf**2 - 23357.D0/135.D0*z2**2*
     &    Ca**2*Cf + 169504.D0/315.D0*z2**3*Cf**3 - 123632.D0/315.D0*
     &    z2**3*Ca*Cf**2 + 7088.D0/63.D0*z2**3*Ca**2*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lqr*as * ( 6*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lqr*as**2 * (  - 34.D0/3.D0*Cf*nf - 
     &    93*Cf**2 + 193.D0/3.D0*Ca*Cf + 48*z3*Cf**2 - 24*z3*Ca*Cf + 16.
     &    D0/3.D0*z2*Cf*nf + 72*z2*Cf**2 - 88.D0/3.D0*z2*Ca*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lqr*as**3 * ( 220.D0/9.D0*Cf*nf**2 + 
     &    230*Cf**2*nf + 1495.D0/2.D0*Cf**3 - 3052.D0/9.D0*Ca*Cf*nf - 
     &    3439.D0/2.D0*Ca*Cf**2 + 3082.D0/3.D0*Ca**2*Cf - 480*z5*Cf**3
     &     + 240*z5*Ca*Cf**2 + 80*z5*Ca**2*Cf - 64.D0/27.D0*z3*Cf*nf**2
     &     - 496.D0/3.D0*z3*Cf**2*nf - 992*z3*Cf**3 + 3440.D0/27.D0*z3*
     &    Ca*Cf*nf + 5368.D0/3.D0*z3*Ca*Cf**2 - 22600.D0/27.D0*z3*Ca**2
     &    *Cf - 608.D0/27.D0*z2*Cf*nf**2 - 272*z2*Cf**2*nf - 720*z2*
     &    Cf**3 + 7504.D0/27.D0*z2*Ca*Cf*nf + 1552*z2*Ca*Cf**2 - 20720.D
     &    0/27.D0*z2*Ca**2*Cf + 704*z2*z3*Cf**3 - 352*z2*z3*Ca*Cf**2 + 
     &    464.D0/5.D0*z2**2*Cf**2*nf + 1968.D0/5.D0*z2**2*Cf**3 - 344.D0
     &    /15.D0*z2**2*Ca*Cf*nf - 2912.D0/5.D0*z2**2*Ca*Cf**2 + 1964.D0/
     &    15.D0*z2**2*Ca**2*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lqr**2*as**2 * ( 2*Cf*nf + 18*Cf**2
     &     - 11*Ca*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lqr**2*as**3 * (  - 68.D0/9.D0*Cf*
     &    nf**2 - 92*Cf**2*nf - 270*Cf**3 + 850.D0/9.D0*Ca*Cf*nf + 551*
     &    Ca*Cf**2 - 2429.D0/9.D0*Ca**2*Cf + 32*z3*Cf**2*nf + 288*z3*
     &    Cf**3 - 16*z3*Ca*Cf*nf - 320*z3*Ca*Cf**2 + 88*z3*Ca**2*Cf + 
     &    32.D0/9.D0*z2*Cf*nf**2 + 48*z2*Cf**2*nf + 144*z2*Cf**3 - 352.D
     &    0/9.D0*z2*Ca*Cf*nf - 264*z2*Ca*Cf**2 + 968.D0/9.D0*z2*Ca**2*
     &    Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lqr**3*as**3 * ( 8.D0/9.D0*Cf*nf**2
     &     + 12*Cf**2*nf + 36*Cf**3 - 88.D0/9.D0*Ca*Cf*nf - 66*Ca*Cf**2
     &     + 242.D0/9.D0*Ca**2*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr*as * (  - 6*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr*as**2 * ( 2.D0/3.D0*Cf*nf + 93*
     &    Cf**2 - 17.D0/3.D0*Ca*Cf - 48*z3*Cf**2 + 24*z3*Ca*Cf + 16.D0/
     &    3.D0*z2*Cf*nf - 72*z2*Cf**2 - 88.D0/3.D0*z2*Ca*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr*as**3 * ( 34.D0/9.D0*Cf*nf**2 - 
     &    275.D0/3.D0*Cf**2*nf - 1495.D0/2.D0*Cf**3 - 40*Ca*Cf*nf + 
     &    2348.D0/3.D0*Ca*Cf**2 + 1657.D0/18.D0*Ca**2*Cf + 480*z5*Cf**3
     &     - 240*z5*Ca*Cf**2 - 80*z5*Ca**2*Cf + 32.D0/9.D0*z3*Cf*nf**2
     &     + 256.D0/3.D0*z3*Cf**2*nf + 992*z3*Cf**3 - 400.D0/9.D0*z3*Ca
     &    *Cf*nf - 4048.D0/3.D0*z3*Ca*Cf**2 + 3104.D0/9.D0*z3*Ca**2*Cf
     &     - 160.D0/27.D0*z2*Cf*nf**2 + 40*z2*Cf**2*nf + 720*z2*Cf**3
     &     + 2672.D0/27.D0*z2*Ca*Cf*nf - 100*z2*Ca*Cf**2 - 8992.D0/27.D0
     &    *z2*Ca**2*Cf - 704*z2*z3*Cf**3 + 352*z2*z3*Ca*Cf**2 + 272.D0/
     &    5.D0*z2**2*Cf**2*nf - 1968.D0/5.D0*z2**2*Cf**3 - 8.D0/5.D0*
     &    z2**2*Ca*Cf*nf - 1136.D0/5.D0*z2**2*Ca*Cf**2 + 4*z2**2*Ca**2*
     &    Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr*lqr*as**2 * (  - 36*Cf**2 )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr*lqr*as**3 * ( 72*Cf**2*nf + 540*
     &    Cf**3 - 420*Ca*Cf**2 - 576*z3*Cf**3 + 288*z3*Ca*Cf**2 - 288*
     &    z2*Cf**3 )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr*lqr**2*as**3 * (  - 12*Cf**2*nf
     &     - 108*Cf**3 + 66*Ca*Cf**2 )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr**2*as**2 * (  - 2*Cf*nf + 18*
     &    Cf**2 + 11*Ca*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr**2*as**3 * ( 4.D0/9.D0*Cf*nf**2
     &     + 20*Cf**2*nf - 270*Cf**3 - 146.D0/9.D0*Ca*Cf*nf - 131*Ca*
     &    Cf**2 + 493.D0/9.D0*Ca**2*Cf - 32*z3*Cf**2*nf + 288*z3*Cf**3
     &     + 16*z3*Ca*Cf*nf + 32*z3*Ca*Cf**2 - 88*z3*Ca**2*Cf + 32.D0/9.
     &    D0*z2*Cf*nf**2 - 48*z2*Cf**2*nf + 144*z2*Cf**3 - 352.D0/9.D0*
     &    z2*Ca*Cf*nf + 264*z2*Ca*Cf**2 + 968.D0/9.D0*z2*Ca**2*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr**2*lqr*as**3 * (  - 12*Cf**2*nf
     &     + 108*Cf**3 + 66*Ca*Cf**2 )
      CHg0tnnnll1 = CHg0tnnnll1 + lfr**3*as**3 * (  - 8.D0/9.D0*Cf*
     &    nf**2 + 12*Cf**2*nf - 36*Cf**3 + 88.D0/9.D0*Ca*Cf*nf - 66*Ca*
     &    Cf**2 - 242.D0/9.D0*Ca**2*Cf )
      CHg0tnnnll1 = CHg0tnnnll1 + 1 

       CHg0tnnnll2=  + as**3 * ( 8*Cf*N4*nfv - 160.D0/3.D0*z5*Cf*N4*nfv
     &     + 28.D0/3.D0*z3*Cf*N4*nfv + 20*z2*Cf*N4*nfv - 4.D0/5.D0*
     &    z2**2*Cf*N4*nfv ) 

       fonnnll1=  + as * (  - 16*Cf + 16*z2*Cf )
      fonnnll1 = fonnnll1 + as**2 * ( 127.D0/6.D0*Cf*nf + 511.D0/4.D0*
     &    Cf**2 - 1535.D0/12.D0*Ca*Cf + 8.D0/9.D0*z3*Cf*nf - 60*z3*
     &    Cf**2 + 604.D0/9.D0*z3*Ca*Cf - 64.D0/3.D0*z2*Cf*nf - 198*z2*
     &    Cf**2 + 376.D0/3.D0*z2*Ca*Cf + 552.D0/5.D0*z2**2*Cf**2 - 92.D0
     &    /5.D0*z2**2*Ca*Cf )
      fonnnll1 = fonnnll1 + as**3 * (  - 7081.D0/243.D0*Cf*nf**2 - 421.D
     &    0/3.D0*Cf**2*nf - 5599.D0/6.D0*Cf**3 + 110651.D0/243.D0*Ca*Cf
     &    *nf + 74321.D0/36.D0*Ca*Cf**2 - 1505881.D0/972.D0*Ca**2*Cf - 
     &    608.D0/9.D0*z5*Cf**2*nf + 1328*z5*Cf**3 - 8*z5*Ca*Cf*nf - 
     &    5512.D0/9.D0*z5*Ca*Cf**2 - 204*z5*Ca**2*Cf + 16.D0/81.D0*z3*
     &    Cf*nf**2 + 9448.D0/27.D0*z3*Cf**2*nf - 460*z3*Cf**3 - 24512.D0
     &    /81.D0*z3*Ca*Cf*nf - 51508.D0/27.D0*z3*Ca*Cf**2 + 139345.D0/
     &    81.D0*z3*Ca**2*Cf + 32*z3**2*Cf**3 + 592.D0/3.D0*z3**2*Ca*
     &    Cf**2 - 400.D0/3.D0*z3**2*Ca**2*Cf + 1072.D0/27.D0*z2*Cf*
     &    nf**2 + 9064.D0/27.D0*z2*Cf**2*nf + 2936.D0/3.D0*z2*Cf**3 - 
     &    44540.D0/81.D0*z2*Ca*Cf*nf - 66544.D0/27.D0*z2*Ca*Cf**2 + 
     &    130295.D0/81.D0*z2*Ca**2*Cf - 256.D0/3.D0*z2*z3*Cf**2*nf - 
     &    400*z2*z3*Cf**3 + 880.D0/9.D0*z2*z3*Ca*Cf*nf + 3680.D0/3.D0*
     &    z2*z3*Ca*Cf**2 - 7228.D0/9.D0*z2*z3*Ca**2*Cf + 448.D0/135.D0*
     &    z2**2*Cf*nf**2 - 36208.D0/135.D0*z2**2*Cf**2*nf - 5972.D0/5.D0
     &    *z2**2*Cf**3 )
      fonnnll1 = fonnnll1 + as**3 * ( 1156.D0/135.D0*z2**2*Ca*Cf*nf + 
     &    258304.D0/135.D0*z2**2*Ca*Cf**2 - 23357.D0/135.D0*z2**2*Ca**2
     &    *Cf + 169504.D0/315.D0*z2**3*Cf**3 - 123632.D0/315.D0*z2**3*
     &    Ca*Cf**2 + 7088.D0/63.D0*z2**3*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lqr*as * ( 6*Cf )
      fonnnll1 = fonnnll1 + lqr*as**2 * (  - 34.D0/3.D0*Cf*nf - 93*
     &    Cf**2 + 193.D0/3.D0*Ca*Cf + 48*z3*Cf**2 - 24*z3*Ca*Cf + 16.D0/
     &    3.D0*z2*Cf*nf + 72*z2*Cf**2 - 88.D0/3.D0*z2*Ca*Cf )
      fonnnll1 = fonnnll1 + lqr*as**3 * ( 220.D0/9.D0*Cf*nf**2 + 230*
     &    Cf**2*nf + 1495.D0/2.D0*Cf**3 - 3052.D0/9.D0*Ca*Cf*nf - 3439.D
     &    0/2.D0*Ca*Cf**2 + 3082.D0/3.D0*Ca**2*Cf - 480*z5*Cf**3 + 240*
     &    z5*Ca*Cf**2 + 80*z5*Ca**2*Cf - 64.D0/27.D0*z3*Cf*nf**2 - 496.D
     &    0/3.D0*z3*Cf**2*nf - 992*z3*Cf**3 + 3440.D0/27.D0*z3*Ca*Cf*nf
     &     + 5368.D0/3.D0*z3*Ca*Cf**2 - 22600.D0/27.D0*z3*Ca**2*Cf - 
     &    608.D0/27.D0*z2*Cf*nf**2 - 272*z2*Cf**2*nf - 720*z2*Cf**3 + 
     &    7504.D0/27.D0*z2*Ca*Cf*nf + 1552*z2*Ca*Cf**2 - 20720.D0/27.D0
     &    *z2*Ca**2*Cf + 704*z2*z3*Cf**3 - 352*z2*z3*Ca*Cf**2 + 464.D0/
     &    5.D0*z2**2*Cf**2*nf + 1968.D0/5.D0*z2**2*Cf**3 - 344.D0/15.D0
     &    *z2**2*Ca*Cf*nf - 2912.D0/5.D0*z2**2*Ca*Cf**2 + 1964.D0/15.D0
     &    *z2**2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lqr**2*as**2 * ( 2*Cf*nf + 18*Cf**2 - 11*Ca
     &    *Cf )
      fonnnll1 = fonnnll1 + lqr**2*as**3 * (  - 68.D0/9.D0*Cf*nf**2 - 
     &    92*Cf**2*nf - 270*Cf**3 + 850.D0/9.D0*Ca*Cf*nf + 551*Ca*Cf**2
     &     - 2429.D0/9.D0*Ca**2*Cf + 32*z3*Cf**2*nf + 288*z3*Cf**3 - 16
     &    *z3*Ca*Cf*nf - 320*z3*Ca*Cf**2 + 88*z3*Ca**2*Cf + 32.D0/9.D0*
     &    z2*Cf*nf**2 + 48*z2*Cf**2*nf + 144*z2*Cf**3 - 352.D0/9.D0*z2*
     &    Ca*Cf*nf - 264*z2*Ca*Cf**2 + 968.D0/9.D0*z2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lqr**3*as**3 * ( 8.D0/9.D0*Cf*nf**2 + 12*
     &    Cf**2*nf + 36*Cf**3 - 88.D0/9.D0*Ca*Cf*nf - 66*Ca*Cf**2 + 242.
     &    D0/9.D0*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lfr*as * (  - 6*Cf )
      fonnnll1 = fonnnll1 + lfr*as**2 * ( 2.D0/3.D0*Cf*nf + 93*Cf**2 - 
     &    17.D0/3.D0*Ca*Cf - 48*z3*Cf**2 + 24*z3*Ca*Cf + 16.D0/3.D0*z2*
     &    Cf*nf - 72*z2*Cf**2 - 88.D0/3.D0*z2*Ca*Cf )
      fonnnll1 = fonnnll1 + lfr*as**3 * ( 34.D0/9.D0*Cf*nf**2 - 275.D0/
     &    3.D0*Cf**2*nf - 1495.D0/2.D0*Cf**3 - 40*Ca*Cf*nf + 2348.D0/3.D
     &    0*Ca*Cf**2 + 1657.D0/18.D0*Ca**2*Cf + 480*z5*Cf**3 - 240*z5*
     &    Ca*Cf**2 - 80*z5*Ca**2*Cf + 32.D0/9.D0*z3*Cf*nf**2 + 256.D0/3.
     &    D0*z3*Cf**2*nf + 992*z3*Cf**3 - 400.D0/9.D0*z3*Ca*Cf*nf - 
     &    4048.D0/3.D0*z3*Ca*Cf**2 + 3104.D0/9.D0*z3*Ca**2*Cf - 160.D0/
     &    27.D0*z2*Cf*nf**2 + 40*z2*Cf**2*nf + 720*z2*Cf**3 + 2672.D0/
     &    27.D0*z2*Ca*Cf*nf - 100*z2*Ca*Cf**2 - 8992.D0/27.D0*z2*Ca**2*
     &    Cf - 704*z2*z3*Cf**3 + 352*z2*z3*Ca*Cf**2 + 272.D0/5.D0*z2**2
     &    *Cf**2*nf - 1968.D0/5.D0*z2**2*Cf**3 - 8.D0/5.D0*z2**2*Ca*Cf*
     &    nf - 1136.D0/5.D0*z2**2*Ca*Cf**2 + 4*z2**2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lfr*lqr*as**2 * (  - 36*Cf**2 )
      fonnnll1 = fonnnll1 + lfr*lqr*as**3 * ( 72*Cf**2*nf + 540*Cf**3
     &     - 420*Ca*Cf**2 - 576*z3*Cf**3 + 288*z3*Ca*Cf**2 - 288*z2*
     &    Cf**3 )
      fonnnll1 = fonnnll1 + lfr*lqr**2*as**3 * (  - 12*Cf**2*nf - 108*
     &    Cf**3 + 66*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lfr**2*as**2 * (  - 2*Cf*nf + 18*Cf**2 + 11
     &    *Ca*Cf )
      fonnnll1 = fonnnll1 + lfr**2*as**3 * ( 4.D0/9.D0*Cf*nf**2 + 20*
     &    Cf**2*nf - 270*Cf**3 - 146.D0/9.D0*Ca*Cf*nf - 131*Ca*Cf**2 + 
     &    493.D0/9.D0*Ca**2*Cf - 32*z3*Cf**2*nf + 288*z3*Cf**3 + 16*z3*
     &    Ca*Cf*nf + 32*z3*Ca*Cf**2 - 88*z3*Ca**2*Cf + 32.D0/9.D0*z2*Cf
     &    *nf**2 - 48*z2*Cf**2*nf + 144*z2*Cf**3 - 352.D0/9.D0*z2*Ca*Cf
     &    *nf + 264*z2*Ca*Cf**2 + 968.D0/9.D0*z2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lfr**2*lqr*as**3 * (  - 12*Cf**2*nf + 108*
     &    Cf**3 + 66*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lfr**3*as**3 * (  - 8.D0/9.D0*Cf*nf**2 + 12
     &    *Cf**2*nf - 36*Cf**3 + 88.D0/9.D0*Ca*Cf*nf - 66*Ca*Cf**2 - 
     &    242.D0/9.D0*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb*as**2 * (  - 224.D0/27.D0*Cf*nf + 1616.
     &    D0/27.D0*Ca*Cf - 56*z3*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb*as**3 * ( 3712.D0/729.D0*Cf*nf**2 + 6*
     &    Cf**2*nf - 125252.D0/729.D0*Ca*Cf*nf - 25856.D0/27.D0*Ca*
     &    Cf**2 + 594058.D0/729.D0*Ca**2*Cf + 384*z5*Ca**2*Cf + 64.D0/9.
     &    D0*z3*Cf*nf**2 + 608.D0/9.D0*z3*Cf**2*nf + 1808.D0/27.D0*z3*
     &    Ca*Cf*nf + 896*z3*Ca*Cf**2 - 24656.D0/27.D0*z3*Ca**2*Cf - 
     &    3584.D0/27.D0*z2*Cf**2*nf + 1648.D0/81.D0*z2*Ca*Cf*nf + 25856.
     &    D0/27.D0*z2*Ca*Cf**2 - 12784.D0/81.D0*z2*Ca**2*Cf - 896*z2*z3
     &    *Ca*Cf**2 + 352.D0/3.D0*z2*z3*Ca**2*Cf + 64.D0/5.D0*z2**2*
     &    Cf**2*nf - 32.D0/5.D0*z2**2*Ca*Cf*nf - 176.D0/5.D0*z2**2*
     &    Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb*lqr*as * (  - 8*Cf )
      fonnnll1 = fonnnll1 + lnNb*lqr*as**2 * ( 80.D0/9.D0*Cf*nf + 128*
     &    Cf**2 - 536.D0/9.D0*Ca*Cf - 128*z2*Cf**2 + 16*z2*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb*lqr*as**3 * (  - 800.D0/81.D0*Cf*nf**2
     &     - 288*Cf**2*nf - 1022*Cf**3 + 16408.D0/81.D0*Ca*Cf*nf + 7006.
     &    D0/3.D0*Ca*Cf**2 - 62012.D0/81.D0*Ca**2*Cf - 640.D0/9.D0*z3*
     &    Cf**2*nf + 480*z3*Cf**3 - 7856.D0/9.D0*z3*Ca*Cf**2 + 352*z3*
     &    Ca**2*Cf + 2816.D0/9.D0*z2*Cf**2*nf + 1584*z2*Cf**3 - 320.D0/
     &    9.D0*z2*Ca*Cf*nf - 19904.D0/9.D0*z2*Ca*Cf**2 + 2144.D0/9.D0*
     &    z2*Ca**2*Cf - 4416.D0/5.D0*z2**2*Cf**3 + 2016.D0/5.D0*z2**2*
     &    Ca*Cf**2 - 352.D0/5.D0*z2**2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb*lqr**2*as**2 * (  - 8.D0/3.D0*Cf*nf - 
     &    48*Cf**2 + 44.D0/3.D0*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb*lqr**2*as**3 * ( 160.D0/27.D0*Cf*nf**2
     &     + 536.D0/3.D0*Cf**2*nf + 744*Cf**3 - 2312.D0/27.D0*Ca*Cf*nf
     &     - 3320.D0/3.D0*Ca*Cf**2 + 7120.D0/27.D0*Ca**2*Cf - 384*z3*
     &    Cf**3 + 192*z3*Ca*Cf**2 - 256.D0/3.D0*z2*Cf**2*nf - 576*z2*
     &    Cf**3 + 32.D0/3.D0*z2*Ca*Cf*nf + 1696.D0/3.D0*z2*Ca*Cf**2 - 
     &    176.D0/3.D0*z2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb*lqr**3*as**3 * (  - 32.D0/27.D0*Cf*
     &    nf**2 - 32*Cf**2*nf - 144*Cf**3 + 352.D0/27.D0*Ca*Cf*nf + 176
     &    *Ca*Cf**2 - 968.D0/27.D0*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb*lfr*as * ( 8*Cf )
      fonnnll1 = fonnnll1 + lnNb*lfr*as**2 * (  - 80.D0/9.D0*Cf*nf - 
     &    128*Cf**2 + 536.D0/9.D0*Ca*Cf + 128*z2*Cf**2 - 16*z2*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb*lfr*as**3 * (  - 32.D0/27.D0*Cf*nf**2
     &     + 288*Cf**2*nf + 1022*Cf**3 - 1672.D0/27.D0*Ca*Cf*nf - 7006.D
     &    0/3.D0*Ca*Cf**2 + 980.D0/3.D0*Ca**2*Cf + 640.D0/9.D0*z3*Cf**2
     &    *nf - 480*z3*Cf**3 - 224.D0/3.D0*z3*Ca*Cf*nf + 7856.D0/9.D0*
     &    z3*Ca*Cf**2 + 176.D0/3.D0*z3*Ca**2*Cf - 2816.D0/9.D0*z2*Cf**2
     &    *nf - 1584*z2*Cf**3 + 320.D0/9.D0*z2*Ca*Cf*nf + 19904.D0/9.D0
     &    *z2*Ca*Cf**2 - 2144.D0/9.D0*z2*Ca**2*Cf + 4416.D0/5.D0*z2**2*
     &    Cf**3 - 2016.D0/5.D0*z2**2*Ca*Cf**2 + 352.D0/5.D0*z2**2*Ca**2
     &    *Cf )
      fonnnll1 = fonnnll1 + lnNb*lfr*lqr*as**2 * ( 96*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb*lfr*lqr*as**3 * (  - 608.D0/3.D0*Cf**2
     &    *nf - 1488*Cf**3 + 3824.D0/3.D0*Ca*Cf**2 + 768*z3*Cf**3 - 384
     &    *z3*Ca*Cf**2 + 1152*z2*Cf**3 - 192*z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb*lfr*lqr**2*as**3 * ( 32*Cf**2*nf + 432
     &    *Cf**3 - 176*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb*lfr**2*as**2 * ( 8.D0/3.D0*Cf*nf - 48*
     &    Cf**2 - 44.D0/3.D0*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb*lfr**2*as**3 * (  - 160.D0/27.D0*Cf*
     &    nf**2 + 24*Cf**2*nf + 744*Cf**3 + 2312.D0/27.D0*Ca*Cf*nf - 
     &    168*Ca*Cf**2 - 7120.D0/27.D0*Ca**2*Cf - 384*z3*Cf**3 + 192*z3
     &    *Ca*Cf**2 + 256.D0/3.D0*z2*Cf**2*nf - 576*z2*Cf**3 - 32.D0/3.D
     &    0*z2*Ca*Cf*nf - 1120.D0/3.D0*z2*Ca*Cf**2 + 176.D0/3.D0*z2*
     &    Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb*lfr**2*lqr*as**3 * ( 32*Cf**2*nf - 432
     &    *Cf**3 - 176*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb*lfr**3*as**3 * ( 32.D0/27.D0*Cf*nf**2
     &     - 32*Cf**2*nf + 144*Cf**3 - 352.D0/27.D0*Ca*Cf*nf + 176*Ca*
     &    Cf**2 + 968.D0/27.D0*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb**2*as * ( 8*Cf )
      fonnnll1 = fonnnll1 + lnNb**2*as**2 * (  - 80.D0/9.D0*Cf*nf - 128
     &    *Cf**2 + 536.D0/9.D0*Ca*Cf + 128*z2*Cf**2 - 16*z2*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb**2*as**3 * ( 800.D0/81.D0*Cf*nf**2 + 
     &    2144.D0/9.D0*Cf**2*nf + 1022*Cf**3 - 16408.D0/81.D0*Ca*Cf*nf
     &     - 17786.D0/9.D0*Ca*Cf**2 + 62012.D0/81.D0*Ca**2*Cf + 640.D0/
     &    9.D0*z3*Cf**2*nf - 480*z3*Cf**3 + 4832.D0/9.D0*z3*Ca*Cf**2 - 
     &    352*z3*Ca**2*Cf - 2816.D0/9.D0*z2*Cf**2*nf - 1584*z2*Cf**3 + 
     &    320.D0/9.D0*z2*Ca*Cf*nf + 19904.D0/9.D0*z2*Ca*Cf**2 - 2144.D0/
     &    9.D0*z2*Ca**2*Cf + 4416.D0/5.D0*z2**2*Cf**3 - 2016.D0/5.D0*
     &    z2**2*Ca*Cf**2 + 352.D0/5.D0*z2**2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb**2*lqr*as**2 * ( 16.D0/3.D0*Cf*nf + 48
     &    *Cf**2 - 88.D0/3.D0*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb**2*lqr*as**3 * (  - 320.D0/27.D0*Cf*
     &    nf**2 - 3968.D0/27.D0*Cf**2*nf - 744*Cf**3 + 4624.D0/27.D0*Ca
     &    *Cf*nf + 23288.D0/27.D0*Ca*Cf**2 - 14240.D0/27.D0*Ca**2*Cf + 
     &    384*z3*Cf**3 + 256*z3*Ca*Cf**2 + 128*z2*Cf**2*nf + 576*z2*
     &    Cf**3 - 64.D0/3.D0*z2*Ca*Cf*nf - 800*z2*Ca*Cf**2 + 352.D0/3.D0
     &    *z2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb**2*lqr**2*as**2 * ( 32*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lqr**2*as**3 * ( 32.D0/9.D0*Cf*
     &    nf**2 - 208.D0/9.D0*Cf**2*nf - 368*Cf**3 - 352.D0/9.D0*Ca*Cf*
     &    nf + 1912.D0/9.D0*Ca*Cf**2 + 968.D0/9.D0*Ca**2*Cf + 512*z2*
     &    Cf**3 - 128*z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lqr**3*as**3 * ( 64.D0/3.D0*Cf**2*
     &    nf + 192*Cf**3 - 352.D0/3.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr*as**2 * (  - 48*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr*as**3 * (  - 208.D0/27.D0*Cf**2
     &    *nf + 744*Cf**3 + 2056.D0/27.D0*Ca*Cf**2 - 384*z3*Cf**3 - 256
     &    *z3*Ca*Cf**2 + 128.D0/3.D0*z2*Cf**2*nf - 576*z2*Cf**3 - 416.D0
     &    /3.D0*z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr*lqr*as**2 * (  - 64*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr*lqr*as**3 * ( 992.D0/9.D0*Cf**2
     &    *nf + 736*Cf**3 - 6992.D0/9.D0*Ca*Cf**2 - 1024*z2*Cf**3 + 256
     &    *z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr*lqr**2*as**3 * (  - 64.D0/3.D0*
     &    Cf**2*nf - 576*Cf**3 + 352.D0/3.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr**2*as**2 * ( 32*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr**2*as**3 * (  - 784.D0/9.D0*
     &    Cf**2*nf - 368*Cf**3 + 5080.D0/9.D0*Ca*Cf**2 + 512*z2*Cf**3
     &     - 128*z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr**2*lqr*as**3 * (  - 64.D0/3.D0*
     &    Cf**2*nf + 576*Cf**3 + 352.D0/3.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**2*lfr**3*as**3 * ( 64.D0/3.D0*Cf**2*
     &    nf - 192*Cf**3 - 352.D0/3.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*as**2 * (  - 32.D0/9.D0*Cf*nf + 176.
     &    D0/9.D0*Ca*Cf )
      fonnnll1 = fonnnll1 + lnNb**3*as**3 * ( 640.D0/81.D0*Cf*nf**2 - 
     &    544.D0/27.D0*Cf**2*nf - 9248.D0/81.D0*Ca*Cf*nf + 4480.D0/27.D0
     &    *Ca*Cf**2 + 28480.D0/81.D0*Ca**2*Cf - 448*z3*Ca*Cf**2 - 512.D0
     &    /9.D0*z2*Cf**2*nf + 128.D0/9.D0*z2*Ca*Cf*nf + 2816.D0/9.D0*z2
     &    *Ca*Cf**2 - 704.D0/9.D0*z2*Ca**2*Cf )
      fonnnll1 = fonnnll1 + lnNb**3*lqr*as**2 * (  - 64*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*lqr*as**3 * (  - 128.D0/27.D0*Cf*
     &    nf**2 + 1088.D0/9.D0*Cf**2*nf + 1024*Cf**3 + 1408.D0/27.D0*Ca
     &    *Cf*nf - 7520.D0/9.D0*Ca*Cf**2 - 3872.D0/27.D0*Ca**2*Cf - 
     &    1024*z2*Cf**3 + 256*z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*lqr**2*as**3 * (  - 64*Cf**2*nf - 
     &    384*Cf**3 + 352*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*lqr**3*as**3 * (  - 256.D0/3.D0*
     &    Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**3*lfr*as**2 * ( 64*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*lfr*as**3 * (  - 1088.D0/9.D0*Cf**2
     &    *nf - 1024*Cf**3 + 7520.D0/9.D0*Ca*Cf**2 + 1024*z2*Cf**3 - 
     &    256*z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*lfr*lqr*as**3 * ( 128.D0/3.D0*Cf**2
     &    *nf + 768*Cf**3 - 704.D0/3.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*lfr*lqr**2*as**3 * ( 256*Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**3*lfr**2*as**3 * ( 64.D0/3.D0*Cf**2*
     &    nf - 384*Cf**3 - 352.D0/3.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**3*lfr**2*lqr*as**3 * (  - 256*Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**3*lfr**3*as**3 * ( 256.D0/3.D0*Cf**3
     &     )
      fonnnll1 = fonnnll1 + lnNb**4*as**2 * ( 32*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**4*as**3 * ( 64.D0/27.D0*Cf*nf**2 - 
     &    640.D0/9.D0*Cf**2*nf - 512*Cf**3 - 704.D0/27.D0*Ca*Cf*nf + 
     &    4288.D0/9.D0*Ca*Cf**2 + 1936.D0/27.D0*Ca**2*Cf + 512*z2*Cf**3
     &     - 128*z2*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**4*lqr*as**3 * ( 640.D0/9.D0*Cf**2*nf
     &     + 192*Cf**3 - 3520.D0/9.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**4*lqr**2*as**3 * ( 256*Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**4*lfr*as**3 * (  - 256.D0/9.D0*Cf**2*
     &    nf - 192*Cf**3 + 1408.D0/9.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**4*lfr*lqr*as**3 * (  - 512*Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**4*lfr**2*as**3 * ( 256*Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**5*as**3 * (  - 256.D0/9.D0*Cf**2*nf
     &     + 1408.D0/9.D0*Ca*Cf**2 )
      fonnnll1 = fonnnll1 + lnNb**5*lqr*as**3 * (  - 256*Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**5*lfr*as**3 * ( 256*Cf**3 )
      fonnnll1 = fonnnll1 + lnNb**6*as**3 * ( 256.D0/3.D0*Cf**3 )
      fonnnll1 = fonnnll1 + 1 

       fonnnll2=  + as**3 * ( 8*Cf*N4*nfv - 160.D0/3.D0*z5*Cf*N4*nfv + 
     &    28.D0/3.D0*z3*Cf*N4*nfv + 20*z2*Cf*N4*nfv - 4.D0/5.D0*z2**2*
     &    Cf*N4*nfv ) 

CC>>>>>>>>>>>>>>>>>>>>>>>>>N3LO
      IF(IMODE .EQ. 1) THEN
        IF(INFV.EQ.0) THEN !! non NFV term
         FUNC = fonnnll1*Zinv                   !!SV
        ELSEIF(INFV.EQ.1) THEN !! NFV term
         FUNC = fonnnll2*Zinv                   !!SV
        ENDIF
      ELSEIF(IMODE .EQ. 2) THEN
         IF (REAL(w).ge. 1.0D0) THEN
          FUNC =(0.0D0,0.0D0)
         ELSEIF(INFV.EQ.0) THEN !! non NFV term
          FUNC = (CHg0tnnnll1*zexp(gnnnll)-fonnnll1)*Zinv   !!RESUM-SV
         ELSEIF(INFV.EQ.1) THEN !! NFV term
          FUNC = (CHg0tnnnll2*zexp(gnnnll)-fonnnll2)*Zinv   !!RESUM-SV
         ENDIF
      ENDIF

        ENDIF
