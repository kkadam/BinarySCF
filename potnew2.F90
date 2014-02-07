      SUBROUTINE SETBDY(NCALL,ISYM)                                             
!C*                                                                              
!C*    THESE ROUTINES WHICH DEAL WITH THE POTENIAL SOLVER                        
!C*    ALLOW FULL 2PI AND ODD MODES (WOODWARD 2/6/89)                            
!C                                                                               
!C...  THIS ROUTINE SIMPLY INITIALIZES VARIOUS ARRAYS BEFORE A CALL TO           
!C     BDYGEN IS MADE.                                                           
!C     IF NCALL = 0, THEN INNITIALIZE ALL ARRAYS IN COMMON BLOCK /BDY/.          
!C              .GT.0, GRID HAS MOVED, SO RE-EVALUATE RBDY, JPOS, AND KPO        
!C                                                                               
!C...  DIMENSION COSM AND SINM (LMAX,10).                                        
!C     THE DIMENSION OF RBDY,JPOS,KPOS, AND IARAY DEPENDS ON SYMMETRIES          
!C     USED (SEE ALSO LAST DIMENSION OF BDYTRM IN ROUTINE BDYGEN):               
!C     IF ABS(ISYM) = 1 OR 8, DIMENSION THEM (2*JMAX + KMAX -1).                 
!C                  = 2,3, OR 9, DIMENSION THEM (JMAX + KMAX -1).                
!C                                                                               
      include "prmtr.h"
!c     PARAMETER (JMAX=64, JMAX1 = JMAX+1, JMAX2 = JMAX+2,                       
!c    &           KMAX=32, KMAX1 = KMAX+1, KMAX2 = KMAX+2,                       
!c    &           LMAX=64, JKM1  = 1*JMAX+KMAX-1)                                
!C                                                                               
      COMMON /BLOK6/   DTHETA,COSIGN(LMAX),SIGN(LMAX),PI,GRAV                   
      COMMON /GRID/    R(JMAX2),Z(KMAX2),RHF(JMAX2),ZHF(KMAX2),        &         
                      G(JMAX2),HH(KMAX2),ROF3N,ZOF3N,A1NEWR,A1NEWZ             
      COMMON /BDY/     COSM(LMAX,10),SINM(LMAX,10),BDYCHK,             &
                      RBDY(JKM1),JPOS(JKM1),KPOS(JKM1),JKMAX            
      DIMENSION        IARAY(JKM1),RB2(JKM1)                                    
!C                                                                               
!C                                                                               
!c     write(6,*)'entering routine setbdy'
!c     write(6,*)ncall,isym
!c     write(6,*)'array zhf'
!c     write(6,1110)r
!c     write(6,1110)z
!c     write(6,1110)rhf
!c     write(6,1110)zhf
      ISYMA = IABS(ISYM)                                                        
      IF(NCALL.NE.0)GO TO 200                                                   
        BDYCHK = 1.0                                                            
        DO 5 I=1,10                                                             
          DO 4 L=1,LMAX                                                         
            COSM(L,I)=0.0                                                       
            SINM(L,I)=0.0                                                       
    4     CONTINUE                                                              
    5   CONTINUE                                                                
        DO 15 L=1,LMAX                                                          
          PSI = FLOAT(L-1)*DTHETA + 0.5*DTHETA                                  
          DO 14 M=1,10                                                          
            PRODCT = PSI*FLOAT(M)                                               
            COSM(L,M)=COS(PRODCT)                                               
            SINM(L,M)=SIN(PRODCT)                                               
   14     CONTINUE                                                              
   15   CONTINUE                                                                
  200 CONTINUE                                                                  
      JKSTOP = JKM1                                                             
!CJWW  JKSTOP=JMAX + KMAX - 1                                                    
!CJWW  IF(ISYMA.EQ.1.OR.ISYMA.EQ.8)JKSTOP=2*JMAX+KMAX-1                          
      DO 205 JK=1,JKSTOP                                                        
        RBDY(JK)=0.0                                                            
        JPOS(JK)=0                                                              
        KPOS(JK)=0                                                              
  205 CONTINUE                                                                  
      I=0                                                                       
      K=KMAX1                                                                   
      ZK2=ZHF(K)**2                                                             
!c
!c
      do 215 J=2,JMAX1
      I=I+1
      RBDY(I)=SQRT(ZK2+RHF(J)**2)
  215 IARAY(I)=I
      IF(ISYMA.NE.1.AND.ISYMA.NE.8)GO TO 228
      K=1
      ZK2=ZHF(1)**2
      DO 225 J=2,JMAX1
      I=I+1
      RBDY(I)=SQRT(ZK2 + RHF(J)**2)
  225 IARAY(I)=I
  228 CONTINUE
      J=JMAX1
      RJ2=RHF(J)**2
      DO 230 K=2,KMAX
      I=I+1
      RBDY(I)=SQRT(RJ2+ZHF(K)**2)
  230 IARAY(I)=I
      JKMAX=I
!C
!C
!c     write(6,*)'jkmax = ',jkmax
!c     write(6,*)'array rbdy before call to sort'
!c     write(6,1110)(rbdy(jk),jk=1,jkmax)
           DO 232 JK=1,JKMAX
  232      RB2(JK)=RBDY(JK)
      CALL SORT(RB2,IARAY,JKMAX)
           DO 234 JK=1,JKMAX
  234      RBDY(JK)=RB2(JK)
!c     write(6,*)'array rbdy after SORT'
!c     write(6,1110)(rbdy(jk),jk=1,jkmax)
 1110 format(1p,8e10.2)
!C  
!c
!C LOOPS 215, 225, 230 MODIFIED FOR VECTOR ... JWW (6/13/89)                     
!c     DO 215 J=2,JMAX1                                                          
!c       RBDY(I+J-1)=SQRT(ZK2+RHF(J)**2)                                         
!c       IARAY(I+J-1)=I+J-1                                                      
!c 215 CONTINUE                                                                  
!c     I=I+JMAX                                                                  
!c     IF(ISYMA.NE.1.AND.ISYMA.NE.8)GO TO 228                                    
!c       K=1                                                                     
!c       ZK2=ZHF(1)**2                                                           
!c       DO 225 J=2,JMAX1                                                        
!c         RBDY(I+J-1)=SQRT(ZK2 + RHF(J)**2)                                     
!c         IARAY(I+J-1)=I+J-1                                                    
!c 225   CONTINUE                                                                
!c       I=I+JMAX                                                                
!c 228 CONTINUE                                                                  
!c     J=JMAX1                                                                   
!c     RJ2=RHF(J)**2                                                             
!c     DO 230 K=2,KMAX                                                           
!c       RBDY(I+K-1)=SQRT(RJ2+ZHF(K)**2)                                         
!c       IARAY(I+K-1)=I+K-1                                                      
!c 230 CONTINUE                                                                  
!c     I=I+KMAX-1                                                                
!c     JKMAX=I                                                                   
!CJWW  I SHOULD EQUAL JKM1 AT THIS POINT                                         
!C                                                                               
!c     CALL SSORTX(RBDY,1,JKMAX,IARAY)                                           
!C                                                                               
!C  ESSL SORT WITH INDEXING FOR 3090 (WOODWARD 1/23/89)                          
!C                                                                               
!c
!c
      JKSKIP = JMAX                                                             
      IF(ISYMA.EQ.1.OR.ISYMA.EQ.8)JKSKIP=2*JMAX                                 
      DO 245 I=1,JKMAX                                                          
        II=IARAY(I)                                                             
        IF(II.GT.JMAX)GO TO 238                                                 
          KPOS(I)=KMAX1                                                         
          JPOS(I)=II+1                                                          
          GO TO 245                                                             
  238   IF(II.GT.JKSKIP)GO TO 240                                               
          KPOS(I)=1                                                             
          JPOS(I)=(II-JMAX) + 1                                                 
          GO TO 245                                                             
  240   JPOS(I)=JMAX1                                                           
        KPOS(I)=(II-JKSKIP) + 1                                                 
  245 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
!C                                                                               
!C ***********************************************************                   
!C                                                                               
      SUBROUTINE BDYGEN(MAXTRM,ISYM,REDGE)                                      
!C                                                                               
!C...   CALL PARAMETERS  ...                                                     
!C  MAXTRM = MAXIMUM L TO BE USED IN SPHERICAL HARMONIC EXPANSION.               
!C           FOR GREATEST EFFICIENCY, USE 4, 8, OR 10.                           
!C           10 IS MAXIMUM ALLOWABLE VALUE.                                      
!C  ISYM = NEGATIVE, ALL MASS TOTALLY WITHIN GRID BOUNDARY R"S.                  
!C       = POSITIVE, USE GENERAL EXPANSION SINCE SOME MASS OUTSIDE.              
!C  ABS(ISYM) = 1,  NO SYMMETRIES.                                               
!C            = 2,  SYM. THRU EQUATORIAL PLANE, FULL 2-PI.                       
!C            = 3,  PI-SYMMETRY AND SYM. THRU EQUATORIAL PLANE.                  
!C            = 8,  2-D WITH NO SYMMETRIES.                                      
!C            = 9,  2-D WITH SYM. THRU EQUATORIAL PLANE.                         
!C  REDGE = 0.0,  MASS CAN BE ANYWHERE IN GRID.                                  
!C        .GT.0.0, MASS ENTIRELY WITHIN SPHERE OF RADIUS = REDGE.                
!C                                                                               
!C                                                                               
      include "prmtr.h"
!c     PARAMETER (JMAX=64, JMAX1 = JMAX+1, JMAX2 = JMAX+2,                       
!c    &           KMAX=32, KMAX1 = KMAX+1, KMAX2 = KMAX+2,                       
!c    &           LMAX=64, JKM1  = 1*JMAX+KMAX-1)                                
!C                                                                               
      COMMON /BLOK6/   DTHETA,COSIGN(LMAX),SIGN(LMAX),PI,GRAV                   
      COMMON /GRID/    R(JMAX2),Z(KMAX2),RHF(JMAX2),ZHF(KMAX2),          &       
                      G(JMAX2),HH(KMAX2),ROF3N,ZOF3N,A1NEWR,A1NEWZ             
      COMMON /POIS/    PHI(JMAX2,KMAX2,LMAX),RHO(JMAX2,KMAX2,LMAX)              
      COMMON /INSIDE/  TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT                      
      COMMON /BDY/     COSM(LMAX,10),SINM(LMAX,10),BDYCHK,               &
                      RBDY(JKM1),JPOS(JKM1),KPOS(JKM1),JKMAX
!C                                                                               
!C                                                                               
      COMMON/TERMS/T00,T10,T11,T20,T21,T22, T30,T31,T32,T33, T40,T41,    &       
      T42,T43,T44, T50,T51,T52,T53,T54,T55, T60,T61,T62,T63,T64,T65,T66  &      
      , T70,T71,T72,T73,T74,T75,T76,T77, T80,T81,T82,T83,T84,T85,T86,    &      
      T87,T88, T90,T91,T92,T93,T94,T95,T96,T97,T98,T99, T100,T101,       &      
      T102,T103,T104,T105,T106,T107,T108,T109,T1010                            
!C                                                                               
!C                                                                               
!C...  THE DIMENSIONS OF TERM,Q,CM,SM,TRMTOT,TRMIN AND TRMOUT                    
!C     SHOULD NEVER BE ALTERED.                                                  
!C     THEY ALLOW FOR EXPANSIONS THROUGH L=10, M=10.                             
!C                                                                               
!C                                                                               
!C...  THE LAST DIMENSION OF ARRAY BDYTRM, HOWEVER, SHOULD BE SET EQUAL          
!C     TO THE SIZE OF ARRAYS RBDY, JPOS, AND KPOS.  (SINCE BDYTRM IS SO          
!C     LARGE, IT CAN, IF NEED BE, BE REDUCED FURTHER IF "MAXTRM.LT.10".          
!C     THE FIRST DIMENSION OF BDYTRM MUST BE .GE."2*NUMTRM", WHERE               
!C     NUMTRM IS DETERMINED IN LOOP 5 OF THIS SUBROUTINE.)                       
!C                                                                               
      DIMENSION TERM(66),Q(10),CM(10),SM(10)                                    
      DIMENSION TRMTOT(132),TRMIN(132),TRMOUT(132),BDYTRM(132,2,JKM1)           
      EQUIVALENCE (TERM(1),T00)                                                 
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
!C      IN A SPHERICAL HARMONIC EXPANSION, TRADITIONALLY,                        
!C                                                                               
!C (1)     YLM(THETA,PSI) = SQRT((2L+1)/4PI*FACTORIAL(L-M)/FACTORIAL(L+M)        
!C                          *PLM(X)*EXP(I*M*PSI),                                
!C                                                                               
!C         YL,-M(THETA,PSI) = (-1)**M*COMPLEX CONJUGATE(YLM),                    
!C                                                                               
!C              WHERE:  X=COS(THETA)   AND   I=SQRT(-1).                         
!C                                                                               
!C      THE EXPANSION OF BOUNDARY POTENTIALS IN TERMS OF SPHERICAL               
!C      HARMONICS LEADS TO PRODUCTS OF YLM AND COMPLEX CONJUGATE(YLM), SO        
!C      A MORE USEFUL DEFINITION OF YLM IS                                       
!C                                                                               
!C (2)     YLM(THETA,PSI) = SQRT((2L+1)/4PI)*SQRT(DLM)*BLM(THETA,PSI)            
!C                          *EXP(I*M*PSI).                                       
!C                                                                               
!C      IN THIS EXPRESSION FOR YLM, PLM(X) HAS BEEN FACTORED INTO A LEADI        
!C      NUMERICAL COEFFICIENT "COEF" AND A THETA-DEPENDENT EXPRESSION BLM        
!C      THEN,                                                                    
!C                                                                               
!C         DLM = (FACTORIAL(L-M)/FACTORIAL(L+M))*COEF**2.                        
!C                                                                               
!C      IN PRACTICE, DLM ALWAYS CONSISTS OF AN ODD INTEGER DIVIDED BY 2**        
!C      THE POWER N VARIES WITH L AND M, BUT FOR TERMS THROUGH L=10, N           
!C      VARIES FROM 0 TO 18.  THE FOLLOWING DATA STATEMENT CONTAINS EXACT        
!C      VALUES OF 2**(-N) FOR N = 5 THROUGH 18; E.G., TW8 = 2**(-8).             
!C      THESE TERMS WILL BE USED TO CALCULATE APPROPRIATE DLM'S.                 
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
      DATA TW5,TW6,TW7,TW8,TW9,TW10,TW11,TW12,TW13,TW14,TW15,TW16,TW17,     &    
      TW18/0.03125,0.015625,7.8125E-3,3.90625E-3,1.953125E-3,               &   
      9.765625E-04,4.8828125E-4,2.44140625E-4,1.220703125E-4,               &   
      6.103515625E-5,3.051757813E-5,1.525878907E-5,7.629394535E-6,          &   
      3.814697266E-6/                                                          
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
!C      THE FOLLOWING STATEMENT FUNCTIONS ARE THE BLM'S USED IN THE              
!C     DEFINITION OF YLM IN EQUATION (2), ABOVE.  VARIABLES THAT ARE            
!C      MULTIPLIED BY NUMERICAL COEFFICIENTS ARE ALWAYS EVEN POWERS OF           
!C      COS(THETA); A VARIABLE PRECEDING A PARENTHETICAL EXPRESSION IS AL        
!C      COS(THETA) TO THE FIRST POWER; A VARIABLE TRAILING A PARENTHETICA        
!C      EXPRESSION IS ALWAYS SOME POWER (EITHER EVEN OR ODD) OF SIN(THETA        
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
      B20(A) = 3.0*A - 1.0                                                      
      B30(A,D) = D*(5.0*A - 3.0)                                                
      B31(A,D) = (5.0*A - 1.0)*D                                                
      B40(A,D) = 35.0*A - 30.0*D + 3.0                                          
      B41(A,D,E) = D*(7.0*A - 3.0)*E                                            
      B42(A,D) = (7.0*A - 1.0)*D                                                
      B50(A,D,E) = E*(63.0*A - 70.0*D + 15.0)                                   
      B51(A,D,E) = (21.0*A - 14.0*D + 1.0)*E                                    
      B52(A,D,E) = D*(3.0*A - 1.0)*E                                            
      B53(A,D) = (9.0*A - 1.0)*D                                                
      B60(A,D,E) = 231.0*A - 315.0*D + 105.0*E - 5.0                            
      B61(A,D,E,F) = E*(33.0*A - 30.0*D + 5.0)*F                                
      B62(A,D,E) = (33.0*A - 18.0*D + 1.0)*E                                    
      B63(A,D,E) = D*(11.0*A - 3.0)*E                                           
      B64(A,D) = (11.0*A - 1.0)*D                                               
      B70(A,D,E,F) = F*(429.0*A - 693.0*D + 315.0*E - 35.0)                     
      B71(A,D,E,F) = (429.0*A - 495.0*D + 135.0*E -5.0)*F                       
      B72(A,D,E,F) = E*(143.0*A - 110.0*D + 15.0)*F                             
      B73(A,D,E) = (143.0*A - 66.0*D + 3.0)*E                                   
      B74(A,D,E) = D*(13.0*A - 3.0)*E                                           
      B75(A,D) = (13.0*A - 1.0)*D                                               
      B80(A,D,E,F) = 6435.0*A - 12012.0*D + 6930.0*E - 1260.0*F + 35.0          
      B81(A,D,E,F,GG) = F*(715.0*A - 1001.0*D +385.0*E - 35.0)*GG               
      B82(A,D,E,F) = (143.0*A - 143.0*D + 33.0*E -1.0)*F                        
      B83(A,D,E,F) = E*(39.0*A - 26.0*D + 3.0)*F                                
      B84(A,D,E) = (65.0*A - 26.0*D + 1.0)*E                                    
      B85(A,D,E) = D*(5.0*A - 1.0)*E                                            
      B86(A,D) = (15.0*A - 1.0)*D                                                
      B90(A,D,E,F,GG) = GG*(12155.0*A - 25740.0*D + 18018.0*E               &               
                      -4620.0*F + 315.0)                                       
      B91(A,D,E,F,GG) = (2431.0*A - 4004.0*D + 2002.0*E - 308.0*F + 7.0)    &       
                      *GG                                                      
      B92(A,D,E,F,GG) = F*(221.0*A - 273.0*D + 91.0*E - 7.0)*GG                 
      B93(A,D,E,F) = (221.0*A - 195.0*D + 39.0*E - 1.0)*F                       
      B94(A,D,E,F) = E*(17.0*A - 10.0*D + 1.0)*F                                
      B95(A,D,E) = (85.0*A - 30.0*D + 1.0)*E                                    
      B96(A,D,E) = D*(17.0*A - 3.0)*E                                           
      B97(A,D) = (17.0*A - 1.0)*D                                               
      B100(A,D,E,F,GG) = 46189.0*A - 109395.0*D + 90090.0*E - 30030.0*F     &    
                       + 3465.0*GG - 63.0                                      
      B101(A,D,E,F,GG,X) = GG*(4199.0*A - 7956.0*D + 4914.0*E - 1092.0*F    &    
                         + 63.0)*X                                             
      B102(A,D,E,F,GG) = (4199.0*A - 6188.0*D + 2730.0*E - 364.0*F          &    
                       + 7.0)*GG                                               
      B103(A,D,E,F,GG) = F*(323.0*A - 357.0*D + 105.0*E - 7.0)*GG               
      B104(A,D,E,F) = (323.0*A - 255.0*D + 45.0*E - 1.0)*F                      
      B105(A,D,E,F) = E*(323.0*A - 170.0*D + 15.0)*F                            
      B106(A,D,E) = (323.0*A - 102.0*D + 3.0)*E                                 
      B107(A,D,E) = D*(19.0*A - 3.0)*E                                          
      B108(A,D) = (19.0*A - 1.0)*D                                              
!C                                                                               
!C                                                                               
!C                                                                               
!C      FINISHED LISTING NEEDED STATEMENT FUNCTIONS.                             
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
!C      THE POTENTIAL AT ANY POINT (RB,THETAB,PSIB) IS GIVEN BY:                 
!C                                                                               
!C         -GRAV*SUM OVER ALL L,M (M.NE.0) OF                                    
!C                                                                               
!C         2.0*COS(M*PSIB)*DLM*BLM(THETAB)                                       
!C            *(RB**-(L+1)*MASSIN1(L,M) + RB**L*MASSOUT1(L,M))                   
!C        +2.0*SIN(M*PSIB)*DLM*BLM(THETAB)                                       
!C            *(RB**-(L+1)*MASSIN2(L,M) + RB**L*MASSOUT2(L,M))                   
!C                                                                               
!C         PLUS -GRAV*SUM OVER ALL L,M=0 OF                                      
!C                                                                               
!C         DL0*BL0(THETAB)*(RB**-(L+1)*MASSIN1(L,0) + RB**L*MASSOUT1(L,0)        
!C                                                                               
!C      GIVEN THAT THE POINT (RB,THETAB,PSIB) IS AT A SPHERICAL RADIUS           
!C      RSPHER, THE TERMS MASSIN AND MASSOUT (EACH A FUNCTION OF L AND M)        
!C      ARE INTEGRALS OVER THE MASS DISTRIBUTION INSIDE AND OUTSIDE,             
!C      RESPECTIVELY, OF RSPHER.   LETTING TERM(R,THETA,PSI) =                   
!C      BLM(THETA)*RHO(R,THETA,PSI)*DVOLUME(R,THETA),                            
!C                                                                               
!C         MASSIN1(L,M) = SUM OVER ALL THETA,PSI,R INSIDE RSPHER OF              
!C                        R**L*COS(M*PSI)*TERM(R,THETA,PSI),                     
!C                                                                               
!C         MASSIN2(L,M) = SUM OVER ALL THETA,PSI,R INSIDE RSPHER OF              
!C                        R**L*SIN(M*PSI)*TERM(R,THETA,PSI),                     
!C                                                                               
!C         MASSOUT1(L,M) = SUM OVER ALL THETA,PSI,R OUTSIDE RSPHER OF            
!C                         R**-(L+1)*COS(M*PSI)*TERM(R,THETA,PSI),               
!C                                                                               
!C         MASSOUT2(L,M) = SUM OVER ALL THETA,PSI,R OUTSIDE RSPHER OF            
!C                         R**-(L+1)*SIN(M*PSI)*TERM(R,THETA,PSI).               
!C                                                                               
!C                                                                               
!C                                                                               
!C                                                                               
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC           
!C                                                                   C           
!C                                                                   C           
!C                        NOW BEGIN PROGRAM.                         C           
!C                                                                   C           
!C                                                                   C           
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC           
!C                                                                               
!C                                                                               
!c     write(6,*)'entering routine bdygen'
      TMASS=0.0                                                                 
      igrid = 0
      IF(BDYCHK.NE.1.0)GO TO 990                                                
      IF(IGRID.GT.0)CALL SETBDY(IGRID,ISYM)                                     
      ISYMA = IABS(ISYM)                                                        
      NTRM = 1                                                                  
      DO 5 LL = 1,MAXTRM                                                        
        NTRM = NTRM + LL + 1                                                    
    5 CONTINUE                                                                  
      NUMTRM = 2*NTRM                                                           
      LELMAX = MAXTRM + 1                                                       
      JMX = JKMAX                                                               
      ISMAX = 2                                                                 
      NUMOUT = NUMTRM                                                           
      IF(ISYM.GT.0) GO TO 6                                                     
        JMX = 1                                                                 
        ISMAX = 1                                                               
        NUMOUT = 1                                                              
    6 CONTINUE                                                                  
      DO 8 I = 1,NUMTRM                                                         
        TRMIN(I) = 0.0                                                          
    8 CONTINUE                                                                  
      DO 9 I = 1,NUMOUT                                                         
        TRMOUT(I) = 0.0                                                         
    9 CONTINUE                                                                  
      DO 12 JK = 1,JMX                                                          
        DO 11 IS = 1,ISMAX                                                      
         DO 10 I = 1,NUMTRM                                                     
             BDYTRM(I,IS,JK) = 0.0                                              
   10     CONTINUE                                                              
   11   CONTINUE                                                                
   12 CONTINUE                                                                  
      F1 = 1.0                                                                  
      IF(ISYMA.EQ.2.OR.ISYMA.EQ.9)F1=0.0                                        
      IF(ISYMA.EQ.3)F1=0.0                                                      
      F2 = 1.0                                                                  
      IF(ISYMA.EQ.2.OR.ISYMA.EQ.9)F2=2.0                                        
      IF(ISYMA.EQ.3)F2=4.0                                                      
!C                                                                               
!C...     BEGIN INTEGRALS OVER MASS DISTRIBUTION.     ...                        
      DO 300 K = 2,KMAX                                                         
        KP = K+1                                                                
        ZZ = ZHF(K)                                                             
        Z2 = ZZ*ZZ                                                              
        ZDEL = Z(KP) - Z(K)                                                     
        DO 300 J = 2,JMAX                                                       
          JP = J+1                                                              
          RR = RHF(J)                                                           
          R2 = RR*RR                                                            
          RSPHER = SQRT(R2 + Z2)                                                
          IF(REDGE.GT.0.0.AND.RSPHER.GT.REDGE)GO TO 300                         
!C                                                                               
!C                                                                               
!C     STEP ONE:  FOR THIS R,THETA, CALCULATE THE PRODUCT                        
!C                BLM(THETA)*R**L FOR ALL L,M.  THIS PRODUCT WILL END UP         
!C                ARRAY TERM, SINCE TERM IS EQUIVALENCED TO T00, T10, ETC        
!C                                                                               
          C = ZZ/RSPHER                                                         
          S = RR/RSPHER                                                         
          C2 = C*C                                                              
          C4 = C2*C2                                                            
          S2 = S*S                                                              
          S3 = S*S2                                                             
          S4 = S2*S2                                                            
          IF(MAXTRM.LE.4)GO TO 602                                              
            C6 = C2*C4                                                          
            C8 = C4*C4                                                          
            S5 = S*S4                                                           
            S6 = S2*S4                                                          
            S7 = S*S6                                                           
            S8 = S4*S4                                                          
            IF(MAXTRM.LE.8)GO TO 602                                            
              C10 = C2*C8                                                       
              S9 = S*S8                                                         
              S10 = S2*S8                                                       
  602     CONTINUE                                                              
!C... FOLLOWING TWO LINES REDUNDANT.  J WOODWARD, JUNE 1988                      
!C         DO 604 LL = 1,MAXTRM                                                  
!C 604     Q(LL) = 0.0                                                           
!C                                                                               
          DO 605 I = 1,NTRM                                                     
            TERM(I) = 0.0                                                       
  605     CONTINUE                                                              
          DO 608 LL = 1,MAXTRM                                                  
            Q(LL) = RSPHER**LL                                                  
  608     CONTINUE                                                              
!C...      EL.EQ.0 THRU 4                                                        
          T00 = 1.0                                                             
          T10 = C*Q(1)*F1                                                       
          T20 = B20(C2)*Q(2)                                                    
          T30 = B30(C2,C)*Q(3)*F1                                               
!C...  *F1 ON T30 COEFF. ADDED BY J WOODWARD, JUNE 1988                          
          T40 = B40(C4,C2) *Q(4)                                                
          IF(ISYMA.GE.8) GO TO 611                                              
            T22 = S2*Q(2)                                                       
            T42 = B42(C2,S2) *Q(4)                                              
            T44 = S4*Q(4)                                                       
            IF(ISYMA.EQ.3) GO TO 611                                            
              T11 = S*Q(1)                                                      
              T31 = B31(C2,S) *Q(3)                                             
              T33 = S3*Q(3)                                                     
              IF(ISYMA.EQ.2) GO TO 611                                          
                T21 = C*S*Q(2)                                                  
                T32 = C*S2*Q(3)                                                 
                T41 = B41(C2,C,S) *Q(4)                                         
                T43 = C*S3*Q(4)                                                 
  611     CONTINUE                                                              
          IF(MAXTRM.LE.4) GO TO 620                                             
!C...      EL.EQ.5 THRU 8                                                        
          T50 = B50(C4,C2,C) *Q(5)*F1                                           
          T60 = B60(C6,C4,C2) *Q(6)                                             
          T70 = B70(C6,C4,C2,C) *Q(7)*F1                                        
          T80 = B80(C8,C6,C4,C2) *Q(8)                                          
          IF(ISYMA.GE.8) GO TO 613                                              
            QH = Q(6)                                                           
            T62 = B62(C4,C2,S2) *QH                                             
            T64 = B64(C2,S4) *QH                                                
            T66 = S6*QH                                                         
            QH = Q(8)                                                           
            T82 = B82(C6,C4,C2,S2) *QH                                          
            T84 = B84(C4,C2,S4) *QH                                             
            T86 = B86(C2,S6) *QH                                                
            T88 = S8*QH                                                         
            IF(ISYMA.EQ.3) GO TO 613                                            
              QH = Q(5)                                                         
              T51 = B51(C4,C2,S) *QH                                            
              T53 = B53(C2,S3) *QH                                              
              T55 = S5*QH                                                       
              QH = Q(7)                                                         
              T71 = B71(C6,C4,C2,S) *QH                                         
              T73 = B73(C4,C2,S3) *QH                                           
              T75 = B75(C2,S5) *QH                                              
              T77 = S7*QH                                                       
              IF(ISYMA.EQ.2) GO TO 613                                          
                T52 = B52(C2,C,S2) *Q(5)                                        
                T54 = C*S4*Q(5)                                                 
                QH = Q(6)                                                       
                T61 = B61(C4,C2,C,S) *QH                                        
                T63 = B63(C2,C,S3) *QH                                          
                T65 = C*S5*QH                                                   
                QH = Q(7)                                                       
                T72 = B72(C4,C2,C,S2) *QH                                       
                T74 = B74(C2,C,S4) *QH                                          
                T76 = C*S6*QH                                                   
                QH = Q(8)                                                       
                T81 = B81(C6,C4,C2,C,S) *QH                                     
                T83 = B83(C4,C2,C,S3) *QH                                       
                T85 = B85(C2,C,S5)*QH                                           
                T87 = C*S7*QH                                                   
  613     CONTINUE                                                              
          IF(MAXTRM.LE.8)GO TO 620                                              
!C...      EL.EQ.9 THRU 10                                                       
          T90 = B90(C8,C6,C4,C2,C) *Q(9)*F1                                     
          T100 = B100(C10,C8,C6,C4,C2) *Q(10)                                   
          IF(ISYMA.GE.8)GO TO 617                                               
            QH = Q(10)                                                          
            T102 = B102(C8,C6,C4,C2,S2) *QH                                     
            T104 = B104(C6,C4,C2,S4) *QH                                        
            T106 = B106(C4,C2,S6) *QH                                           
            T108 = B108(C2,S8) *QH                                              
            T1010 = S10*QH                                                      
            IF(ISYMA.EQ.3) GO TO 617                                            
              QH = Q(9)                                                         
              T91 = B91(C8,C6,C4,C2,S)*QH                                       
              T93 = B93(C6,C4,C2,S3)*QH                                         
              T95 = B95(C4,C2,S5) *QH                                           
              T97 = B97(C2,S7) *QH                                              
              T99 = S9*QH                                                       
              IF(ISYMA.EQ.2) GO TO 617                                          
                T92 = B92(C6,C4,C2,C,S2) *QH                                    
                T94 = B94(C4,C2,C,S4) *QH                                       
                T96 = B96(C2,C,S6) *QH                                          
                T98 = C*S8*QH                                                   
                QH = Q(10)                                                      
                T101 = B101(C8,C6,C4,C2,C,S) *QH                                
                T103 = B103(C6,C4,C2,C,S3) *QH                                  
                T105 = B105(C4,C2,C,S5) *QH                                     
                T107 = B107(C2,C,S7) *QH                                        
                T109 = C*S9*QH                                                  
  617     CONTINUE                                                              
  620     CONTINUE                                                              
!C                                                                               
!C                                                                               
!C                                                                               
!C      STEP TWO:  CALCULATE MASSIN1 AND MASSIN2.  THESE MASSIN TERMS WIL        
!C                 DEPEND ON L AND M, AND IN GENERAL ON THE SPHERICAL RAD        
!C                 OF THE BOUNDARY ZONE BEING CONSIDERED.  FOR EACH L,M          
!C                 TERM (I=1,NUMTRM,2), AND FOR EACH BOUNDARY CELL               
!C                 (JK=1,JKMAX), MASSIN1 IS STORED IN BDYTRM(I,1,JK)             
!C                 AND MASSIN2 IS STORED IN BDYTRM(I+1,1,JK).                    
!C                                                                               
!C                                                                               
!C...       DEPENDING ON CHOSEN SYMMETRY, F2 = 1.0,2.0 OR 4.0 TO ACCOUNT         
!C          FOR ALL MASS.                                                        
          DV = 0.5*DTHETA*ZDEL*(R(JP)**2 - R(J)**2)*F2                          
!c         write(6,*)'dtheta,zdel,f2 = ',dtheta,zdel,f2
!c         write(6,*)'dv = ',dv	
          DO 204 M = 1,MAXTRM                                                   
            CM(M) = 0.0                                                         
            SM(M) = 0.0                                                         
  204     CONTINUE                                                              
          SUMASS = 0.0                                                          
          DO 211 L = 1,LMAX                                                     
!C...      H = MASS IN SINGLE GRID CELL (DV HAS MASS SYMMETRIES BUILT IN         
            H = RHO(J,K,L)*DV                                                   
!C         if(j.eq.k)then
!C         write(6,*)'rho = ',rho(j,k,l)
!C         endif
            SUMASS = SUMASS + H                                                 
            DO 210 M = 1,MAXTRM                                                 
              CM(M) = CM(M) + COSM(L,M)*H                                       
              SM(M) = SM(M) + SINM(L,M)*H                                       
  210       CONTINUE                                                            
  211     CONTINUE                                                              
!c         if(j.eq.k)then
!c         write(6,*)'sumass = ',sumass
!c         write(6,*)'values of sumass along diagonal inside bdygen'
!c         write(6,1111)sumass
!c         endif
!CJWW                 DO 213 I = 1,NUMTRM                                        
!CJWW         213     TRMTOT(I) = 0.0                                            
          NCNT = 0                                                              
          DO 225 LEL = 1,LELMAX                                                 
            LL = LEL - 1                                                        
            MMAX = LEL                                                          
            DO 224 MM = 1,MMAX                                                  
              M = MM -1                                                         
              NCNT = NCNT + 1                                                   
              IP = 2*NCNT                                                       
              I = IP -1                                                         
              IF(M.EQ.0) GO TO 216                                              
                TRMTOT(I) = TERM(NCNT)*CM(M)                                    
                TRMTOT(IP) = TERM(NCNT)*SM(M)                                   
                GO TO 224                                                       
  216         CONTINUE                                                          
              TRMTOT(I) = TERM(NCNT)*SUMASS                                     
  224       CONTINUE                                                            
  225     CONTINUE                                                              
          TMASS=TMASS+TRMTOT(1)
!C  STORE SUM OF TERMS, TO BE USED LATER FOR EACH BNDRY CELL.                    
          IF(ISYM.LT.0) GO TO 1250                                              
!C...               RSPHER**L HAS BEEN USED IN EXPANSIONS.                       
!C*****             VERY DIFFERENT STUFF INSERTED.                               
          JKSTRT= 1                                                             
          IF(RSPHER.LE.RBDY(1))GO TO 227                                        
            JKSSS = 1                                                           
            DO 233 JK=2,JKMAX                                                   
              IF(RSPHER.LE.RBDY(JK))GO TO 234                                   
              JKSSS = JK                                                        
  233       CONTINUE                                                            
            WRITE(6,181)J,K,RSPHER,RBDY(JKMAX),JKMAX                            
  181       FORMAT(//,' STOPSTOPSTOPSPOPSTOPSTOP',//,                     &               
           ' INSIDE BDYGEN, AFTER LOOP 233,',/,                           &     
           ' RSPHER(',I3,',',I3,') =',1P,E13.5,', WHICH IS GREATER',/,    &     
           ' THAN RBDY(JKMAX) =',1P,E13.5,', SO STOP THE PROGRAM....',    &     
           /,' NOTE:  JKMAX =',I5,//,' STOPSTOPSTOPSTOPSTOPSTOP',//)           
            STOP                                                                
  234       JKSTRT = JKSSS +1                                                   
  227     CONTINUE                                                              
!C                                                                               
!C...        FOR JK=JKSTRT,JKMAX MASS RING IS INTERIOR TO RBDY(JK),              
!C           SO STORE RESULT IN RBDY(JKMAX)'S BDYTRM FOR USE BY ALL              
!C           BOUNDARY CELLS IN LOOP 312, BELOW.                                  
!C           CHANGE LOOP TO VECTOR ADDITION (WOODWARD 1/23/89)                   
!C                                                                               
          DO 228 I=1,NUMTRM                                                     
            BDYTRM(I,1,JKMAX) = BDYTRM(I,1,JKMAX)+TRMTOT(I)                     
228       CONTINUE                                                              
          IF(JKSTRT .EQ. 1)GO TO 300                                            
!C                                                                               
!C... FOR JK=1,JKSSS MASS RING IS EXTERIOR TO RBDY(JK), SO                       
!C      STORE RESULT IN RBDY(I,1,JK) SO THAT IT CAN BE SUBTRACTED                
!C      FROM RBDY(JKMAX)'S BDYTRM IN LOOP 312, BELOW.                            
          DO 239 JK=1,JKSSS                                                     
            DO 235 I=1,NUMTRM                                                   
              BDYTRM(I,1,JK) = BDYTRM(I,1,JK) + TRMTOT(I)                       
 235        CONTINUE                                                            
 239      CONTINUE                                                              
!C                                                                               
!C...       TERMS FOR SAME RING OF MASS SHOULD BE DONE AGAIN IF ITS MASS         
!C          LIES OUTSIDE SOME BOUNDARY RADII.                                    
!C                                                                               
!C                                                                               
!C      STEP THREE:  CALCULATE MASSOUT1 AND MASSOUT2.  FOR EACH L,M TERM         
!C      (OPTIONAL)   (I=1,NUMTRM,2) AND FOR EACH BOUNDARY CELL (JK=1,JKMA        
!C                   MASSOUT1 IS STORED IN BDYTRM(I,2,JK) AND MASSOUT2 IS        
!C                   STORED IN BDYTRM(I+1,2,JK).                                 
!C                      NOTE:  FOR ANY PARTICULAR RING OF MASS AT SPHERIC        
!C                   RADIUS R, ITS CONTRIBUTION TO MASSOUT IS JUST               
!C                   R**-(2L+1) TIMES ITS CONTRIBUTION TO MASSIN.                
!C                                                                               
!C                                                                               
          TRMTOT(1) = TRMTOT(1)/RSPHER                                          
          NCNT = 2                                                              
          DO 240 LL = 1,MAXTRM                                                  
            IPOWR = -(2*LL + 1)                                                 
            MMAX = 2*(LL + 1)                                                   
            QH = RSPHER**IPOWR                                                  
            DO 240 MM = 1,MMAX                                                  
              NCNT = NCNT + 1                                                   
              TRMTOT(NCNT) = TRMTOT(NCNT)*QH                                    
  240     CONTINUE                                                              
!C...       RSPHER**-(L+1) HAS BEEN USED IN EXPANSIONS.                          
!C...       MASS RING IS OUTSIDE RBDY(JK).                                       
          DO 246 JK = 1,JKSSS                                                   
           DO 245 I = 1,NUMTRM                                                  
              BDYTRM(I,2,JK) = BDYTRM(I,2,JK) + TRMTOT(I)                       
 245       CONTINUE                                                             
 246      CONTINUE                                                              
          GO TO 300                                                             
 1250     CONTINUE                                                              
!C...       ALL MASS INSIDE GRID BNDRY R"S, SO BDYTRM SAME FOR ALL CELLS.        
          DO 1255 I = 1,NUMTRM                                                  
            BDYTRM(I,1,1) = BDYTRM(I,1,1) + TRMTOT(I)                           
 1255     CONTINUE                                                              
  300   CONTINUE                                                                
  301 CONTINUE                                                                  
!C                                                                               
!C                                                                               
!C      FINISHED CALCULATING MASSIN1, MASSIN2, MASSOUT1, AND MASSOUT2.           
!C                                                                               
!C                                                                               
!C      NOW, FOR EACH BOUNDARY CELL, CALCULATE SPECIFIC VALUE OF THE POTE        
!C                                                                               
!C                                                                               
!c     write(6,*)'maxtrm,isym,redge,tmass = ',maxtrm,isym,redge,tmass
!c     write(6,1111)(rho(jj1,jj1,1),jj1=1,jmax2)
 1111 format(1p,8e10.2)
      DO 400 JK = 1,JKMAX                                                       
        J = JPOS(JK)                                                            
        K = KPOS(JK)                                                            
        RSPHER = RBDY(JK)                                                       
!C                                                                               
!C                                                                               
!C      STEP FOUR:  CALCULATE THE PRODUCT DLM*BLM AND STORE RESULTS IN           
!C                  ARRAY TERM.                                                  
!C                                                                               
!C                                                                               
        C = ZHF(K)/RSPHER                                                       
        S = RHF(J)/RSPHER                                                       
        C2 = C*C                                                                
        C4 = C2*C2                                                              
        S2 = S*S                                                                
        S3 = S*S2                                                               
        S4 = S2*S2                                                              
        IF(MAXTRM.LE.4)GO TO 304                                                
          C6 = C2*C4                                                            
          C8 = C4*C4                                                            
          S5 = S*S4                                                             
          S6 = S2*S4                                                            
          S7 = S*S6                                                             
          S8 = S4*S4                                                            
          IF(MAXTRM.LE.8)GO TO 304                                              
            C10 = C2*C8                                                         
            S9 = S*S8                                                           
            S10 = S2*S8                                                         
  304   CONTINUE                                                                
        DO 305 I = 1,NTRM                                                       
          TERM(I) = 0.0                                                         
  305   CONTINUE                                                                
!C... EL.EQ.0 THRU 4                                                             
        T00 = 1.0                                                               
        T10 = C*F1                                                              
        T20 = 0.25*B20(C2)                                                      
        T30 = 0.25*B30(C2,C)*F1                                                 
        T40 = TW6 *B40(C4,C2)                                                   
        IF(ISYMA.GE.8)GO TO 307                                                 
          T22 = 0.375*S2                                                        
          T42 = 5.0*TW5 *B42(C2,S2)                                             
          T44 = 35.0*TW7 *S4                                                    
          IF(ISYMA.EQ.3)GO TO 307                                               
            T11 = 0.5*S                                                         
            T31 = 0.1875*B31(C2,S)                                              
            T33 = 0.3125*S3                                                     
            IF(ISYMA.EQ.2) GO TO 307                                            
              T21 = 1.5*C*S                                                     
              T32 = 1.875*C*S2                                                  
              T41 = 0.3125 *B41(C2,C,S)                                         
              T43 = 2.1875 *C*S3                                                
  307   CONTINUE                                                                
        IF(MAXTRM.LE.4)GO TO 310                                                
!C...       EL.EQ.5 THRU 8                                                       
        T50 = TW6 *B50(C4,C2,C) *F1                                             
        T60 = TW8 *B60(C6,C4,C2)                                                
        T70 = TW8 *B70(C6,C4,C2,C) *F1                                          
        T80 = TW14 *B80(C8,C6,C4,C2)                                            
        IF(ISYMA.GE.8)GO TO 308                                                 
          T62 = 105.0*TW10 *B62(C4,C2,S2)                                       
          T64 = 63.0*TW9 *B64(C2,S4)                                            
          T66 = 231.0*TW10 *S6                                                  
          T82 = 315.0*TW12 *B82(C6,C4,C2,S2)                                    
          T84 = 693.0*TW13 *B84(C4,C2,S4)                                       
          T86 = 429.0*TW12 *B86(C2,S6)                                          
          T88 = 6435.0*TW15 *S8                                                 
          IF(ISYMA.EQ.3)GO TO 308                                               
            T51 = 15.0*TW7 *B51(C4,C2,S)                                        
            T53 = 35.0*TW8 *B53(C2,S3)                                          
            T55 = 63.0*TW8 *S5                                                  
            T71 = 7.0*TW11 *B71(C6,C4,C2,S)                                     
            T73 = 21.0*TW11 *B73(C4,C2,S3)                                      
            T75 = 231.0*TW11 *B75(C2,S5)                                        
            T77 = 429.0*TW11 *S7                                                
            IF(ISYMA.EQ.2)GO TO 308                                             
              T52 = 105.0*TW5*B52(C2,C,S2)                                      
              T54 = 315.0*TW7 *C*S4                                             
              T61 = 21.0*TW7 *B61(C4,C2,C,S)                                    
              T63 = 105.0*TW8 *B63(C2,C,S3)                                     
              T65 = 693.0*TW8*C*S5                                              
              T72 = 21.0*TW10 *B72(C4,C2,C,S2)                                  
              T74 = 231.0*TW9 *B74(C2,C,S4)                                     
              T76 = 3003.0*TW10 *C*S6                                           
              T81 = 9.0*TW11 *B81(C6,C4,C2,C,S)                                 
              T83 = 1155.0*TW11 *B83(C4,C2,C,S3)                                
              T85 = 9009.0*TW11 *B85(C2,C,S5)                                   
              T87 = 6435.0*TW11 *C*S7                                           
  308   CONTINUE                                                                
        IF(MAXTRM.LE.8)GO TO 310                                                
!C...       EL.EQ.9 THRU 10                                                      
        T90 = TW14 *B90(C8,C6,C4,C2,C)*F1                                       
        T100 = TW16 *B100(C10,C8,C6,C4,C2)                                      
        IF(ISYMA.GE.8)GO TO 309                                                 
          T102 = 165.0*TW17 *B102(C8,C6,C4,C2,S2)                               
          T104 = 2145.0*TW15 *B104(C6,C4,C2,S4)                                 
          T106 = 2145.0*TW18 *B106(C4,C2,S6)                                    
          T108 = 12155.0*TW17 *B108(C2,S8)                                      
          T1010 = 46189.0*TW18 *S10                                             
          IF(ISYMA.EQ.3) GO TO 309                                              
            T91 = 45.0*TW15 *B91(C8,C6,C4,C2,S)                                 
            T93 = 1155.0*TW14 *B93(C6,C4,C2,S3)                                 
            T95 = 1287.0*TW14 *B95(C4,C2,S5)                                    
            T97 = 6435.0*TW16 *B97(C2,S7)                                       
            T99 = 12155.0*TW16*S9                                               
            IF(ISYMA.EQ.2)GO TO 309                                             
              T92 = 495.0*TW12 *B92(C6,C4,C2,C,S2)                              
              T94 = 45045.0*TW13 *B94(C4,C2,C,S4)                               
              T96 = 2145.0*TW12 *B96(C2,C,S6)                                   
              T98 = 109395.0*TW15 *C*S8                                         
              T101 = 55.0*TW15 *B101(C8,C6,C4,C2,C,S)                           
              T103 = 2145.0*TW14 *B103(C6,C4,C2,C,S3)                           
              T105 = 429.0*TW14 *B105(C4,C2,C,S5)                               
              T107 = 36465.0*TW16 *B107(C2,C,S7)                                
              T109 = 230945.0*TW16 *C*S9                                        
  309   CONTINUE                                                                
  310   CONTINUE                                                                
!C                                                                               
!C                                                                               
!C      STEP FIVE:  COMBINE MASSIN (CALLED TRMIN HERE) AND MASSOUT (CALLE        
!C                  TRMOUT) TERMS APPROPRIATELY, KEEPING COSINE DEPENDENT        
!C                  TERMS (TRMTOT(I=ODD)) SEPARATE FROM SINE DEPENDENT TE        
!C                  (TRMTOT(I=EVEN)).                                            
!C                                                                               
!C                                                                               
        IF(ISYM.LT.0)GO TO 1321                                                 
        DO 312 I = 1,NUMTRM                                                     
         TRMIN(I) = BDYTRM(I,1,JKMAX)-BDYTRM(I,1,JK)                            
        TRMOUT(I) = BDYTRM(I,2,JK)                                              
  312   CONTINUE                                                                
        NCNT = 0                                                                
        DO 320 LEL = 1,LELMAX                                                   
          LL = LEL - 1                                                          
          MMAX = LEL                                                            
          IF(LL.EQ.0)GO TO 314                                                  
            RIN = 1.0/RSPHER**MMAX                                              
            ROUT = RSPHER**LL                                                   
            GO TO 315                                                           
  314     CONTINUE                                                              
          RIN = 1.0/RSPHER                                                      
          ROUT = 1.0                                                            
  315     CONTINUE                                                              
          DO 319 MM = 1,MMAX                                                    
            NCNT = NCNT + 1                                                     
            IP = NCNT*2                                                         
            I = IP - 1                                                          
            B = TERM(NCNT)*RIN                                                  
            H = TERM(NCNT)*ROUT                                                 
            TRMTOT(I) = B*TRMIN(I) + H*TRMOUT(I)                                
            TRMTOT(IP) = B*TRMIN(IP) + H*TRMOUT(IP)                             
  319     CONTINUE                                                              
  320   CONTINUE                                                                
        GO TO 331                                                               
 1321   CONTINUE                                                                
!C...       ASSUMING NO MASS OUTSIDE GRID RADIUS.                                
        DO 1323 I = 1,NUMTRM                                                    
          TRMIN(I) = BDYTRM(I,1,1)                                              
 1323   CONTINUE                                                                
        NCNT = 0                                                                
        DO 1330 LEL = 1,LELMAX                                                  
          LL = LEL - 1                                                          
          MMAX = LEL                                                            
          RIN = 1.0/RSPHER**MMAX                                                
          DO 1329 MM = 1,MMAX                                                   
            NCNT = NCNT + 1                                                     
            IP = 2*NCNT                                                         
            I = IP -1                                                           
            B = TERM(NCNT)*RIN                                                  
            TRMTOT(I) = B*TRMIN(I)                                              
            TRMTOT(IP) = B*TRMIN(IP)                                            
 1329     CONTINUE                                                              
 1330   CONTINUE                                                                
  331   CONTINUE                                                                
!C                                                                               
!C                                                                               
!C      FINALLY, STEP SIX:  SUM OVER ALL L,M TERMS, TAKING INTO ACCOUNT          
!C                 THE PSIB ANGLE DEPENDENCE.                                    
!C                                                                               
!C                                                                               
        DO 350 L = 1,LMAX                                                       
          DO 335 M = 1,MAXTRM                                                   
            CM(M) = COSM(L,M)                                                   
            SM(M) = SINM(L,M)                                                   
  335     CONTINUE                                                              
          NCNT = 0                                                              
          SUM = 0.0                                                             
          DO 345 LEL = 1,LELMAX                                                 
            LL = LEL - 1                                                        
            MMAX = LEL                                                          
            DO 345 MM = 1,MMAX                                                  
              M = MM - 1                                                        
              NCNT = NCNT + 2                                                   
              I = NCNT -1                                                       
              IP = NCNT                                                         
              IF(M.EQ.0)GO TO 340                                               
                SUM = SUM + 2.0*(CM(M)*TRMTOT(I) + SM(M)*TRMTOT(IP))            
                GO TO 341                                                       
  340         CONTINUE                                                          
              SUM = SUM + TRMTOT(I)                                             
  341        CONTINUE                                                           
  345       CONTINUE                                                            
  346     CONTINUE                                                              
          PHI(J,K,L) = - GRAV*SUM                                               
  350   CONTINUE                                                                
  400 CONTINUE                                                                  
!C                                                                               
!C                                                                               
!C      FINISHED HARD WORK.  NOW TIDY UP BOUNDARY SYMMETRIES.                    
!C                                                                               
!C                                                                               
      IF(ISYMA.EQ.1)GO TO 550                                                   
!C...  FOR 2-D OR PI-SYM PROBLEMS, SET BOUNDARY CONDITIONS ON Z-AXIS.            
      IF(ISYMA.NE.3.AND.ISYMA.NE.8.AND.ISYMA.NE.9)GO TO 525                     
      DO 505 L=1,LMAX                                                           
  505 PHI(1,KMAX1,L) = PHI(2,KMAX1,L)                                           
      IF(ISYMA.NE.8)GO TO 525                                                   
      DO 510 L=1,LMAX                                                           
  510 PHI(1,1,L) = PHI(2,1,L)                                                   
      GO TO 580                                                                 
!C                                                                               
!C...  FOR SYM THROUGH EQUATORIAL PLANE, SET BNDRY CONDITION THRU PLANE.         
  525 IF(ISYMA.NE.2.AND.ISYMA.NE.3.AND.ISYMA.NE.9)GO TO 550                     
      DO 530 L=1,LMAX                                                           
  530 PHI(JMAX1,1,L) = PHI(JMAX1,2,L)                                           
      IF(ISYMA.NE.2)GO TO 580                                                   
!C                                                                               
!C...  IF NO SYM THRU Z-AXIS, EQUATE CORRECT PHI'S AT J=1.                       
  550 LHAF = LMAX/2                                                             
      DO 560 L=1,LHAF                                                           
        LP=L+LHAF                                                               
        PHI(1,KMAX1,L) = PHI(2,KMAX1,LP)                                        
        PHI(1,KMAX1,LP) = PHI(2,KMAX1,L)                                        
        IF(ISYMA.NE.1)GO TO 560                                                 
        PHI(1,1,L) = PHI(2,1,LP)                                                
        PHI(1,1,LP) = PHI(2,1,L)                                                
  560 CONTINUE                                                                  
  580 CONTINUE                                                                  
      RETURN                                                                    
!C                                                                               
!C                                                                               
!C                                                                               

	!!wth??
	
!  990 WRITE(6,101)BDYCHK                                                        
!  101 FORMAT(//,' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',        
!     1      5X,'BDYCHK =',1P,E10.2,/,5X,'SETBDY HAS NOT BEEN CALLED.',/,        
!     2      5X,'IT MUST BE CALLED AT LEAST ONCE IN ORDER TO INITIALIZE',        
!     3      5X,'THE ARRAYS IN COMMON BLOCK /BDY/',//,                           
!     4      ' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',//)   
	
	
  990 WRITE(6,101)BDYCHK                                                        
  101 FORMAT(//,' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',     &   
           5X,'BDYCHK =',1P,E10.2,/,5X,'SETBDY HAS NOT BEEN CALLED.',/,      &  
           5X,'IT MUST BE CALLED AT LEAST ONCE IN ORDER TO INITIALIZE',      &  
           5X,'THE ARRAYS IN COMMON BLOCK /BDY/',//,                         &  
           ' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',//)  	
	
	
      STOP                                                                      
      END                                                                       
!C *******************************************************************           
      SUBROUTINE POT3(NPOINT,IPRINT,ISYM)                                       
!C                                                                               
!C...  THE CALL PARAMETER ISYM IS DEFINED IN ROUTINE BDYGEN.                     
!C     THIS SUBROUTINE DOES NOT PRESENTLY ALLOW FOR ABS(ISYM)=1 OR 2             
!C                                                                               
!C !!! CHANGES MADE TO ALLOW ISYMA=1,2 BY JOHN WOODWARD 8/10/88 !!!              
!C                                                                               
!c.....lower case modifications introduced by Joel in August 93
!c     to produce a serial poisson solver for Kim Barker's WD merger
!c     project.
      include "prmtr.h"
!c     PARAMETER (JMAX=64, JMAX1 = JMAX+1, JMAX2 = JMAX+2, IJJ = JMAX-1,         
!c    &           KMAX=32, KMAX1 = KMAX+1, KMAX2 = KMAX+2, IKK = KMAX-1,         
!c    &           LMAX=64, LMAX2 = 2*LMAX, L2P1 = LMAX/2 + 1,                    
!c    &           MWFW=650, KKK = KMAX/2 + 1, 
!c    &           ifft=lmax)                                   
!C            (TO DETERMINE SIZE OF MWFW:  LET K = INT(LOG2(IKK)) + 1;           
!C             SET L = 2**(K+1);  THEN MWFW.GE.((K-2)*L+K+5+MAX(2N,6M)).)        
!C            (SEE SUBROUTINE SETUP FOR VALUE OF MWSAVE.)                        
!C                                                                               
      COMMON /BLOK6/   DTHETA,COSIGN(LMAX),SIGN(LMAX),PI,GRAV                   
      COMMON /GRID/    R(JMAX2),Z(KMAX2),RHF(JMAX2),ZHF(KMAX2),           &       
                      G(JMAX2),H(KMAX2),ROF3N,ZOF3N,A1NEWR,A1NEWZ              
      COMMON /POIS/    PHI(JMAX2,KMAX2,LMAX),RHO(JMAX2,KMAX2,LMAX)              
      COMMON /COEFS/   COEF(JMAX1,LMAX,2)                                       
!C     COMMON /SETPHI/  WSAVE(MWSAVE)                                            
!C            (ARRAY A1 IS USED IN FFT;                                          
!C             DIMENSION IT 2*LMAX FOR PI-SYM RUN.)                              
      DIMENSION        A1(ifft),a2(ifft),wsave(2*ifft + 15)                     
!c     COMPLEX*8        B1(LMAX+1)                                               
!C  B1 USED IN FFT (WOODWARD 1/23/89)                                            
!C            (THESE NEXT ARRAYS ARE USED IN BLKTRI.)                            
      DIMENSION        AN(IKK),BN(IKK),CN(IKK),                           &      
                      AM(IJJ),BM(IJJ),CM(IJJ),Y(IJJ,IKK),WFW(MWFW)             
!C            (THESE NEXT ARRAYS STORE QUANTITIES USED TO FIND BLKTRI            
!C             COEFFICIENTS.  ARRAY C STORES LAPLACIAN'S ANGULAR                 
!C             OPERATOR.)                                                        
      DIMENSION        C(IJJ),DENOMR(JMAX),RD2(JMAX),RD3(JMAX),           &      
                      DENOMZ(KMAX),ZD2(KMAX),XLAMM(L2P1)                       
!C  ADD AUX VARS FOR FFT.                                                        
!c     REAL*8 AUX1,AUX2,AUX3                                                     
!c     COMMON/AUX/AUX1(1000),AUX2(1000),AUX3(1000)                               
      DATA ZERO/0.0E0/                                                          
!C                                                                               
!C                                                                               
!C                                                                               
!c     write(6,*)'entering routine pot3'
!c     write(6,*)npoint,iprint,isym
!c     write(6,*)jmax,kmax,lmax
!c     write(6,*)jmax2,kmax2,ifft
!c     write(6,1110)(phi(jjj,kmax1,1),jjj=1,jmax2)
!c     write(6,1110)(phi(jmax1,kk2,1),kk2=1,kmax2)
 1110 format(1p,8e10.2)
      ISYMA = IABS(ISYM)                                                        
      KEQ = 2                                                                   
      IF (ISYMA .EQ. 1 .OR. ISYMA .EQ. 8) KEQ = KKK                             
      N = 2*LMAX                                                                
      IF(ISYMA.EQ.1 .OR. ISYMA .EQ. 2) N=LMAX                                   
!C  ###         THIS 'IF' IS TO ALLOW ODD MODE CALCULATIONS (NON-PI SYM)         
!C  ###         WOODWARD 8/10/88                                                 
      NL2 = LMAX/2                                                              
      IF(LMAX.EQ.1) N=1                                                         
      N1 = NL2 + 1                                                              
!CJWW           SF1=0.5                                                          
!CJWW           SF2=1.0/FLOAT(N)                                                 
      SF1=1.0                                                                   
      SF2=1.0/float(n)                                                        
!c     SCALE = 1.0/FLOAT(N)                                                      
!C              SCALE FACTORS = 1.0 FOR ESSL FFT (WOODWARD 1/23/89)              
!C...      ZERO ARRAY A1.                                                        
      DO 20 L = 1,N                                                             
        A1(L)=0.0                                                               
   20 CONTINUE                                                                  
!C...  PUT,FOR CONVENIENCE, ALL RHO'S INTO PHI ARRAY -- LEAVE BNDRY PHI'S        
!C     ALONE.                                                                    
      DO 30 L=1,LMAX                                                            
        DO 29 K=2,KMAX                                                          
           DO 29 J=2,JMAX                                                       
             PHI(J,K,L)=RHO(J,K,L)                                              
   29   CONTINUE                                                                
   30 CONTINUE                                                                  
      IF(LMAX.EQ.1) GO TO 200                                                   
!C...  NOW, FOR ALL J,K, CALCULATE FOURIER TRANSFORM OF RHO'S AND BNDRY          
!C     PHI'S.                                                                    
!C INITIALIZE THE ESSL FFT (WOODWARD 1/23/89)                                    
!c     CALL SRCFT(1,A1,NDUMX,B1,NDUMY,N,1,1,1.0,                                 
!c    &           AUX1,1000,AUX2,1000,AUX3,1000)                                 
!c     write(6,*)'calling rffti first time'
      call rffti(n,wsave)
      DO 40 K=1,KMAX1                                                           
        DO 40 J=2,JMAX1                                                         
!C...  PUT DENSITIES INTO ARRAY A1, ONE RING AT A TIME.                          
          DO 35 L=1,LMAX                                                        
            A1(L)=PHI(J,K,L)                                                    
            if(isyma.eq.3)A1(L+LMAX) = PHI(J,K,L)                               
   35     CONTINUE                                                              
!CJWW THE A1(L+LMAX)... IS REALLY FOR EVEN MODES ONLY, BUT NO HARM               
!CJWW DONE LEAVING IT IN ALWAYS.                                                 
!C...  PERFORM FOURIER TRANSFORMATION                                            
!CJWW  CALL RFFT(A1,N,1)                                                         
!C  USE ESSL FFT                                                                 
!c       CALL SRCFT(0,A1,NDUMX,B1,NDUMY,N,1,1,1.0,                               
!c    &             AUX1,1000,AUX2,1000,AUX3,1000)                               
!c     write(6,*)'calling rfftf'     
        call rfftf(n,a1,wsave)
        do 737 i=1,n
          a2(i) = a1(i)
  737   continue  
!c       DO 36 I=1,N/2                                                           
!c         I1 = 2*I - 1                                                          
!c         I2 = 2*I                                                              
!c         A1(I1)=REAL(B1(I))                                                    
!c         A1(I2)=AIMAG(B1(I))                                                   
        a1(2) = a2(n)
        do 738 i=3,n
          a1(i) = a2(i-1)
  738   continue  
36      CONTINUE                                                                
!C     FOR IBM3090 FFT ROUTINE, THE COSINE COEFFICIENT OF THE                    
!C     HIGHEST ORDER MODE IS RETURNED AS REAL(B1(LMAX+1)).                       
!C     TO MAKE EVERYTHING CONSISTENT WITH PREVIOUS TREATMENT,                    
!C     PUT THIS RESULT IN A1(2)  (TOHLINE)                                       
!c       A1(2) = REAL(B1(N/2+1))                                                 
!C                                                                               
!C                                                                               
!C...  NOW PUT COSINE COEFFICIENTS INTO PHI(L=1,N1) AND SINE                     
!C     COEFFICIENTS INTO PHI(L=N1+1,LMAX).   REMEMBER FOR PI-SYM THAT            
!C     ONLY EVEN MODES ARE NON-ZERO.                                             
        PHI(J,K,1)=A1(1)*SF1                                                    
        PHI(J,K,N1)=A1(2)*SF1                                                   
        IF (ISYMA .EQ. 1 .OR. ISYMA .EQ. 2) THEN                                
          DO 371 I=2,NL2                                                        
            PHI(J,K,I) = A1(2*I-1)                                              
            PHI(J,K,I+NL2) = A1(2*I)                                            
 371      CONTINUE                                                              
        ELSE                                                                    
          DO 37 I=2,NL2                                                         
            PHI(J,K,I)=A1(4*I-3)*SF1                                            
            PHI(J,K,I+NL2)=A1(4*I-2)*SF1                                        
   37     CONTINUE                                                              
        ENDIF                                                                   
      IF(K.NE.KEQ)GO TO 1045                                                    
!C...       FROM A1 ARRAY, FIGURE OUT AMPLITUDE AND PHASES OF VARIOUS            
!C          MODES.  FOR PI-SYM, ONLY EVEN MODES ARE NON-ZERO.                    
        X1=ABS(A1(1))*SF1                                                       
        IF(X1.EQ.0.0)GO TO 1045                                                 
        MODE=1                                                                  
        COEF(J,MODE,1)=X1                                                       
        COEF(J,MODE,2)=0.0                                                      
        NSTOP=N-1                                                               
        NSTSTP=4                                                                
        IF (ISYMA.EQ.1 .OR. ISYMA.EQ.2) NSTSTP=2                                
!C ### ABOVE 2 LINES ADDED TO ALLOW ODD MODES... VALUES ARE USED IN              
!C ### FOLLOWING LOOP (1034)   (WOODWARD 8/10/88)                                
!C TO ALLOW FOR VECTOR, MAKE PHILAG AN ARRAY OF DIM (LMAX2)                      
        JJJ = NSTOP / NSTSTP                                                    
        DO 1034 III=1,JJJ                                                       
!C         MODE=MODE+1                                                           
!C III + 1 = MODE #                                                              
!C          AC IS COSINE COEFFICIENT                                             
          AC=A1(III*NSTSTP+1)                                                   
!C          AS IS SINE COEFFICIENT                                               
          AS=A1(III*NSTSTP+2)                                                   
          AMP=2.0*SQRT(AC**2 + AS**2)/X1                                        
          IF(AC.NE.0.0)PHILAG=ATAN(-AS/AC)                                      
          IF(AC.EQ.0.0)PHILAG=PI*00.5                                           
          IF(AC.EQ.0.0.AND.AS.EQ.0.0)PHILAG=4.0*PI                              
          IF(AC.LT.0.0)PHILAG=PHILAG + PI                                       
          IF(PHILAG.LT.0.0)PHILAG=PHILAG + 2.0*PI                               
          COEF(J,III+1,1)=AMP                                                   
          COEF(J,III+1,2)=PHILAG*180./PI                                        
 1034   CONTINUE                                                                
        MLAST=LMAX                                                              
        COEF(J,MLAST,1)=2.0*A1(2)/X1                                            
        COEF(J,MLAST,2)=0.0                                                     
 1045   CONTINUE                                                                
!C...  FINISHED AMPLITUDES AND PHASE ANGLES.                                     
   40 CONTINUE                                                                  
!c
!c
!c
  200 CONTINUE                                                                  
!C...  THROUGH WITH FOURIER TRANSFORMATIONS;  NOW SOLVE 2-D EQUATIONS.           
!C... SET UP VARIOUS OPERATORS.                                                  
      DTHET2=1.0/(DTHETA*DTHETA)                                                
      PIG4=4.0*PI*GRAV                                                          
      RD2(1)=0.5/RHF(2)                                                         
      RD3(1)=1.0/RHF(1)                                                         
      DO 210 J=2,JMAX                                                           
        DENOMR(J)=1.0/(RHF(J+1)-RHF(J-1))                                       
        RD2(J)=1.0/(RHF(J+1)-RHF(J))                                            
        RD3(J)=1.0/RHF(J)                                                       
  210 CONTINUE                                                                  
!C     ZD2(1)=0.5/ZHF(2)                                                         
      ZD2(1) = 1.0/(ZHF(2) - ZHF(1))                                            
      DO 220 K=2,KMAX                                                           
        DENOMZ(K)=1.0/(ZHF(K+1)-ZHF(K-1))                                       
        ZD2(K)=1.0/(ZHF(K+1)-ZHF(K))                                            
  220 CONTINUE                                                                  
      DO 230 J=2,JMAX                                                           
        I=J-1                                                                   
        CM(I)=2.*(0.5*RD3(J)+RD2(J))*DENOMR(J)                                  
        AM(I)=2.*(RD2(J-1)-0.5*RD3(J))*DENOMR(J)                                
        C(I)=RD3(J)*RD3(J)*DTHET2                                               
  230 CONTINUE                                                                  
      DO 240 K=2,KMAX                                                           
        I=K-1                                                                   
        CN(I)=2.*ZD2(K)*DENOMZ(K)                                               
        AN(I)=2.*ZD2(K-1)*DENOMZ(K)                                             
        BN(I)=-CN(I)-AN(I)                                                      
  240 CONTINUE                                                                  
!C                                                                               
!C                                                                               
!C                                                                               
!C...  THE BOUNDARY CONDITIONS ON THE POISSON EQUATION ARE ALL HIDDEN            
!C          IN HOW YOU SPECIFY THE FIRST AND LAST VALUES OF ARRAYS               
!C          CM, AM, BM, CN, AN, AND BN.                                          
!C     DIRICHLET CONDITION AT JMAX.....                                          
      IMAX=JMAX-1                                                               
      CMMAX=CM(IMAX)                                                            
      CM(IMAX)=0.0                                                              
!C          CMMAX IS USED IN LOOP 330 AND AT STATEMENT 331 BELOW TO              
!C          COMPLETE CONDITIION.                                                 
!C                                                                               
!C     DIRICHLET CONDITION AT KMAX.....                                          
      IMAX=KMAX-1                                                               
      CNMAX=CN(IMAX)                                                            
      CN(IMAX)=0.0                                                              
!C          CNMAX IS USED IN LOOP 320 BELOW TO COMPLETE CONDITION.               
!C                                                                               
!C     NEUMANN CONDITION AT Z-AXIS (J=2).....                                    
      AMMIN=AM(1)                                                               
      AM(1)=0.0                                                                 
!C          SEE MODIFICATION TO BM(1) AT STATEMENT #341 BELOW.                   
!C                                                                               
      IF(ISYMA.NE.1.AND.ISYMA.NE.8)GO TO 242                                    
!C     DIRICHLET CONDITION AT K=2.....                                           
      ANMIN=AN(1)                                                               
      AN(1)=0.0                                                                 
!C          ANMIN IS USED IN LOOP 325 BELOW TO COMPLETE CONDITION.               
!C                                                                               
      IF(ISYMA.EQ.1.OR.ISYMA.EQ.8)GO TO 245                                     
  242 CONTINUE                                                                  
!C     NEUMANN CONDITION AT K=2 (SYM. THRU EQUATORIAL PLANE).....                
      BN(1)=BN(1)+AN(1)                                                         
      AN(1)=0.0                                                                 
  245 CONTINUE                                                                  
!C                                                                               
!C                                                                               
!C                                                                               
      DO 250 M1=1,N1                                                            
        M = 2*(M1-1)                                                            
        IF (ISYMA.EQ.1 .OR. ISYMA.EQ.2) M=M1-1                                  
!C ###            ABOVE LINE ALLOWS ODD MODES (WOODWARD 8/10/88)                 
        EM=FLOAT(M)                                                             
        XLAMM(M1)=COS(EM*DTHETA)                                                
  250 CONTINUE                                                                  
!C...  COEFFICIENTS BM AND Y WILL VARY WITH M, SO THEY HAVEN'T BEEN              
!C     CALCULATED YET.                                                           
!C...  NOW FOR EACH VALUE OF L, THRU LMAX, CALCULATE PHI'S IN TRANSFORMED        
!C     SPACE.                                                                    
      DO 500 L=1,LMAX                                                           
        IF(L.LE.N1) M1=L                                                        
        IF(L.GT.N1) M1=L + 1 - N1                                               
        POWRM = (-1)**(M1-1)                                                    
        DO 311 K=2,KMAX                                                         
          IK=K-1                                                                
          DO 310 J=2,JMAX                                                       
            IJ=J-1                                                              
!C...  REMEMBER, DENSITIES ARE IN ARRAY PHI RIGHT NOW.                           
            Y(IJ,IK)=PIG4*PHI(J,K,L)                                            
  310     CONTINUE                                                              
  311   CONTINUE                                                                
!C...  NOW CALCULATE BM'S.                                                       
        JSTOP = JMAX-2                                                          
        DO 340 J=2,JSTOP                                                        
          BM(J)=-CM(J)-AM(J)+2.0*(XLAMM(M1)-1.0)*C(J)                           
  340   CONTINUE                                                                
!C...  MODIFY Y AND BM ACCORDING TO THE CHOSEN BOUNDARY CONDITIONS               
!C                                                                               
!C     DIRICHLET CONDITION AT JMAX.....                                          
        J=JMAX                                                                  
        J1=J+1                                                                  
        IJ=J-1                                                                  
        DO 330 K=2,KMAX                                                         
          IK=K-1                                                                
          Y(IJ,IK) = Y(IJ,IK) - CMMAX*PHI(J1,K,L)                               
  330   CONTINUE                                                                
        I=JMAX-1                                                                
  331   BM(I)=-CMMAX-AM(I)+2.0*(XLAMM(M1)-1.0)*C(I)                             
!C                                                                               
!C     DIRICHLET CONDITION AT KMAX.....                                          
        K=KMAX                                                                  
        K1=K+1                                                                  
        IK=K-1                                                                  
        DO 320 J=2,JMAX                                                         
          IJ=J-1                                                                
          Y(IJ,IK) = Y(IJ,IK) - CNMAX*PHI(J,K1,L)                               
  320   CONTINUE                                                                
!C                                                                               
!C     NEUMANN CONDITION AT Z-AXIS (J=2).....                                    
  341   BM(1) = -CM(1) + 2.0*(XLAMM(M1)-1.0)*C(1)                               
        IF (ISYMA.EQ.1 .OR. ISYMA.EQ.2) BM(1)=BM(1)+(POWRM-1.0)*AMMIN           
!C ### ABOVE LINE ALLOWS ODD MODES... (WOODWARD 8/10/88)                         
!C                                                                               
        IF(ISYMA.NE.1.AND.ISYMA.NE.8)GO TO 326                                  
!C     DIRICHLET CONDITION AT K=2.....                                           
        IK=1                                                                    
        DO 325 J=2,JMAX                                                         
          IJ=J-1                                                                
          Y(IJ,IK) = Y(IJ,IK) - ANMIN*PHI(J,1,L)                                
  325   CONTINUE                                                                
  326   CONTINUE                                                                
!C                                                                               
!C                                                                               
!C                                                                               
!C...  THIS COMPLETES THE SET UP OF COEFFICIENTS FOR GIVEN M.                    
        IJ=JMAX-1                                                               
        IK=KMAX-1                                                               
        IF(L.EQ.1) THEN                                                         
!C     write(6,*)'calling blktri(0)'
          CALL BLKTRI(0,1,IK,AN,BN,CN,1,IJ,AM,BM,CM, IJ,Y,             &         
                     IERROR,WFW)                                               
          IF(IERROR.NE.0) WRITE(6,2001)IERROR                                   
        ENDIF                                                                   
 2001   FORMAT(5X,'IERROR = ',I10)                                              
!C     write(6,*)'calling blktri'
        CALL BLKTRI(1,1,IK,AN,BN,CN,1,IJ,AM,BM,CM,IJ,Y,IERROR,WFW)              
        IF(IERROR.NE.0) WRITE(6,2001)IERROR                                     
!C...  SOLUTION ON 2-D GRID IS COMPLETE.                                         
!c
!C...  PUT TRANSFORMED PHI'S FROM Y INTO PHI ARRAY.                              
        DO 351 K=2,KMAX                                                         
          IK=K-1                                                                
          DO 350 J=2,JMAX                                                       
            IJ=J-1                                                              
            PHI(J,K,L)=Y(IJ,IK)                                                 
  350     CONTINUE                                                              
  351   CONTINUE                                                                
!C...  NOW DO ANOTHER VALUE OF M.                                                
  500 CONTINUE                                                                  
!c
!c
!c
      IF(LMAX.EQ.1) GO TO 550                                                   
!C...  ALL TRANSFORMED PHI'S HAVE BEEN CALCULATED.  NOW OBTAIN REAL              
!C     PHI'S BY FOURIER ANALYSIS.                                                
!C  FIRST INITIALIZE FFT (WOODWARD 1/23/89)                                      
!c     write(6,*)'calling rffti second time'
      call rffti(n,wsave)
!c     CALL SCRFT(1,B1,NDUMX,A1,NDUMY,N,1,-1,SCALE,                              
!c    &           AUX1,1000,AUX2,1000,AUX3,1000)                                 
      DO 541 K=1,KMAX1                                                          
        DO 540 J=2,JMAX1                                                        
          DO 505 L = 1,N                                                        
            A1(L)=0.0                                                           
  505     CONTINUE                                                              
          A1(1)=PHI(J,K,1)                                                      
          A1(2)=PHI(J,K,N1)                                                     
!C FIRST SETTLE ODD MODE PROBLEM (FIRST LOOP) THEN EVEN (SECOND)                 
          IF(ISYMA.EQ.1 .OR. ISYMA.EQ.2) THEN                                   
            DO 510 I=2,NL2                                                      
              A1(2*I-1) = PHI(J,K,I)                                            
              A1(2*I)   = PHI(J,K,I+NL2)                                        
  510       CONTINUE                                                            
          ELSE                                                                  
            DO 5100 I=2,NL2                                                     
              A1(4*I-3)=PHI(J,K,I)                                              
              A1(4*I-2)=PHI(J,K,I+NL2)                                          
 5100       CONTINUE                                                            
          ENDIF                                                                 
!C...                                                                            
          do 837 i=1,n
            a2(i) = a1(i)
  837     continue
          do 838 i=3,n
            a1(i-1) = a2(i)
  838     continue
            a1(n) = a2(2)
!c     write(6,*)'calling rfftb'
          call rfftb(n,a1,wsave)
!c         CALL SCRFT(0,B1,NDUMX,A1,NDUMY,N,1,-1,SCALE,                          
!c    &               AUX1,1000,AUX2,1000,AUX3,1000)                             
!C...                                                                            
!C                                                                               
          DO 520 L=1,LMAX                                                       
            PHI(J,K,L)=A1(L)*SF2                                                
  520     CONTINUE                                                              
  540   CONTINUE                                                                
  541 CONTINUE                                                                  
  550 CONTINUE                                                                  
!C                                                                               
!C                                                                               
!C...  CALL ZAXPHI TO DETERMINE PHI ON Z-AXIS.                                   
      CALL ZAXPHI(10,0,ISYM)                                                    
!C                                                                               
!C                                                                               
!C...  CALCULATION OF PHI'S IS NOW FINISHED.                                     
!C ### CHANGES HERE ALLOW ODD MODES (WOODWARD 8/10/88)                           
!C ### CHANGES ARE : IF, LOOPS 554,552,CHANGE IN FIRST 570 L LOOP                
!C ### (USED TO BE DO 580....), ADDITION OF 575 L LOOP.                          
!cjc  Bug that always did pi symmetry on the boundary
!cjc      IF (ISYMA.NE.1 .OR. ISYMA.NE.2) GO TO 562                                 
      IF (ISYMA.NE.1 .OR. ISYMA.NE.2) GO TO 562                                 
      DO 554 L=1,NL2                                                            
        DO 552 K=1,KMAX1                                                        
!cjc  Bug that stepped in L when it should have stepped in K!!!!
!cjc          PHI(1,K,L)=PHI(2,L,L+NL2)                                             
          PHI(1,K,L)=PHI(2,K,L+NL2)                                             
          PHI(1,K,L+NL2)=PHI(2,K,L)                                             
  552   CONTINUE                                                                
  554 CONTINUE                                                                  
      GOTO 572                                                                  
  562 CONTINUE                                                                  
      DO 571 L=1,LMAX                                                           
        DO 570 K=1,KMAX1                                                        
          PHI(1,K,L)=PHI(2,K,L)                                                 
  570   CONTINUE                                                                
  571 CONTINUE                                                                  
  572 CONTINUE                                                                  
      IF(ISYMA.EQ.1.OR.ISYMA.EQ.8)GO TO 576                                     
      DO 575 L=1,LMAX                                                           
        DO 574 J=1,JMAX2                                                        
          PHI(J,1,L)=PHI(J,2,L)                                                 
  574   CONTINUE                                                                
  575 CONTINUE                                                                  
  576 CONTINUE                                                                  
  580 CONTINUE                                                                  
      RETURN                                                                    
!C 990 WRITE(6,101)ISYM                                                          
!C 101 FORMAT(//,' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',        
!C    1      5X,'ISYM =',I5,/,'POT3 DOES NOT YET ALLOW FULL 2PI CALCULATI        
!C    2ON.',//,                                                                  
!C    4      ' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',//)         
!C     RETURN                                                                    
      END                                                                       
!C *************************************************************                 
      SUBROUTINE ZAXPHI(NPOINT,IPRINT,ISYM)                                     
!C                                                                               
      include "prmtr.h"
      parameter (ij = jmax-1, ik = kmax-1, l2 = lmax)
!c     PARAMETER (JMAX=64, JMAX1 = JMAX+1, JMAX2 = JMAX+2, IJ = JMAX-1,          
!c    1           KMAX=32, KMAX1 = KMAX+1, KMAX2 = KMAX+2, IK = KMAX-1,          
!c    2           LMAX=64, L2 = LMAX)                                            
!C     (IT LOOKS LIKE L2 ONLY NEEDS TO = LMAX/2 IF DOING FULL 2*PI,              
!C      BUT IT SHOULDN'T HURT TO LEAVE IT AS LMAX ALWAYS.)                       
!C                                                                               
      COMMON /BLOK6/   DTHETA,COSIGN(LMAX),SIGN(LMAX),PI,GRAV                   
      COMMON /GRID/    R(JMAX2),Z(KMAX2),RHF(JMAX2),ZHF(KMAX2),          &       
                      G(JMAX2),H(KMAX2),ROF3N,ZOF3N,A1NEWR,A1NEWZ              
      COMMON /POIS/    PHI(JMAX2,KMAX2,LMAX),RHO(JMAX2,KMAX2,LMAX)              
      COMMON /COEFS/   COEF(JMAX1,LMAX,2)                                       
      COMMON /INSIDE/  TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT                      
      DIMENSION        X(10),PHI2(10),PH(L2),PHAC(KMAX2)                        
!C            (DIMENSION OF X AND PHI2 SET BY NUMBER OF POINTS IN                
!C             INTERPOLATION SCHEME.......HERE N=10.)                            
!C                                                                               
!C                                                                               
  100 FORMAT(///,' YOU ARE LIMITED TO A 10-POINT INTERPOLATION HERE IN SUBROUTINE ZAXPHI.' &
  ,/,' BUT YOU PUT NPOINT =',I4,///)                      
  101 FORMAT(1P,11E11.3)                                                        
  102 FORMAT(20X,I5,' COEFS=',1P,8E12.3,/,25X,' PHASE=',0P,8F12.2)              
!C     WITH X AND PHI2 DIMENSIONED 10, 10-POINT INTERPOLATION IS THE             
!C     MAXIMUM YOU CAN DO.                                                       
!c     write(6,*)'entering zaxphi'
      ISYMA=IABS(ISYM)                                                          
      IF(NPOINT.GT.10) WRITE(6,100) NPOINT                                      
      IF(NPOINT.GT.10) NPOINT=10                                                
      XPOINT=0.0                                                                
      LM=LMAX                                                                   
      N=NPOINT                                                                  
      ICHK=0                                                                    
      IF(MOD(N,2).NE.0) ICHK=1                                                  
      N=N-ICHK                                                                  
      NUM=N/2                                                                   
      ISPECL=NUM+NUM+1                                                          
      JSP=NUM+1                                                                 
      DO 10 J=2,JSP                                                             
      I1=NUM+(J-1)                                                              
      I2=NUM-(J-2)                                                              
      X(I1)=RHF(J)                                                              
      X(I2)=-RHF(J)                                                             
   10 CONTINUE                                                                  
      IF(ICHK.EQ.1) X(ISPECL)=RHF(NUM+2)                                        
      DO 50 K=2,KMAX1                                                           
      PHMAX=0.0                                                                 
!C...  FOR GIVEN K, CALCULATE PHI ON Z-AXIS AT ALL ANGLES.  STORE                
!C     RESULTS IN PH(L) ARRAY AND IN PHI(JMAX2,K,L).                             
      LHAF=LMAX/2                                                               
      DO 20 L=1,LM                                                              
      LP=L                                                                      
      IF(ISYMA.EQ.1.OR.ISYMA.EQ.2)LP=L+LHAF                                     
      IF(LP.GT.LMAX)LP=L-LHAF                                                   
      DO 15 J=2,JSP                                                             
      I1=NUM+(J-1)                                                              
      I2=NUM-(J-2)                                                              
      PHI2(I1)=PHI(J,K,L)                                                       
      PHI2(I2)=PHI(J,K,LP)                                                      
   15 CONTINUE                                                                  
      IF(ICHK.EQ.1) PHI2(ISPECL)=PHI(NUM+2,K,L)                                 
      IM=NPOINT-1                                                               
      DO 18 I=1,IM                                                              
      P1=PHI2(I)                                                                
      XINV=1.0/X(I)                                                             
      XR=XPOINT*XINV                                                            
      IST=I+1                                                                   
      DO 18 J=IST,NPOINT                                                        
      XRATIO=X(J)*XINV                                                          
      P2=PHI2(J)                                                                
      PHI2(J)=(P1*(XRATIO-XR)+P2*(XR-1.0))/(XRATIO-1.0)                         
   18 CONTINUE                                                                  
      PH(L)=PHI2(NPOINT)                                                        
      IF(ABS(PH(L)).GT.ABS(PHMAX))PHMAX=PH(L)                                   
      PHI(JMAX2,K,L)=PH(L)                                                      
      PHI(JMAX2,K,LP)=PH(L)                                                     
   20 CONTINUE                                                                  
!C...  AT THIS K, FIND MAXIMUM DEVIATION IN PH(L)'S; PUT RESULT IN               
!C     PHAC(K).                                                                  
      ERR=0.0                                                                   
      DO 30 L=1,LM                                                              
      ER=1.0-PH(L)/PHMAX                                                        
      IF(ABS(ER).GT.ABS(ERR))ERR=ER                                             
   30 CONTINUE                                                                  
      PHAC(K)=ERR                                                               
   50 CONTINUE                                                                  
!C...  THE NEXT LOOP FINDS LARGEST OF ALL DEVIATIONS IN PHI'S ON Z-AXIS,         
!C     PUTS VALUE IN PHICHK AND K-VALUE IN KLOCAT.                               
      ERR=0.0                                                                   
      DO 60 K=2,KMAX1                                                           
      XX=ABS(PHAC(K))                                                           
      IF(XX.LE.ABS(ERR)) GO TO 60                                               
      ERR=PHAC(K)                                                               
      KLOCAT=K                                                                  
   60 CONTINUE                                                                  
      PHICHK=ERR                                                                
      IF(IPRINT.NE.0) WRITE(6,101)(PHAC(K),K=2,KMAX1)                           
      IF(LMAX.EQ.1)RETURN                                                       
      IF(IPRINT.NE.0) WRITE(6,102)(J,((COEF(J,M,III),M=1,8),III=1,2),J=12,JMAX)
      
      RETURN                                                                    
      END                                                                       
      FUNCTION EPMACH(DUM)
      EPMACH = 100.0/2.**64.0
      RETURN
      END
      SUBROUTINE SORT(KEY,INDX,NUM)
!C...   DOUBLE SHELL SORT FROM H.M.MURPHY 1967
!C...  J.E.TOHLINE GOT THIS FROM S.H.HODSON 10/17/80.
      REAL KEY(1)
      INTEGER INDX(1),TI
   10 IF(NUM.LT.2)RETURN
      I=1        
   20 I=I+I      
      IF(I.LE.NUM)GO TO 20
      M=I-1      
   30 M=M/2      
      IF(M.LT.1)RETURN
      K=NUM-M    
      DO 50 J=1,K
      I=J        
   40 IM=I+M     
      IF(KEY(I).LE.KEY(IM))GO TO 50
      TK=KEY(I)  
      TI=INDX(I) 
      KEY(I)=KEY(IM)
      INDX(I)=INDX(IM)
      KEY(IM)=TK 
      INDX(IM)=TI
      I=I-M      
      IF(I.GE.1)GO TO 40
   50 CONTINUE   
      GO TO 30   
      END        

