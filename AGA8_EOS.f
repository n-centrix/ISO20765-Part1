     
!ISO 20765 Part 1 - Coded by R.Djigouadi
!www.n-centrix.com
!Release Version Oct-2020  
!Reference:
!BS EN ISO 20765-1:2005
!Natural gas — Calculation of thermodynamic properties Part 1: Gas phase properties for transmission and distribution applications
!
!BS EN ISO 12213-2:2009
!Natural gas — Calculation of compression factor Part 2: Calculation using molarcomposition analysis (ISO 12213-2:2006)
!
!Added: 
!Calculation of isentropic temperature exponent as a derivated properties. 
 
        include 'AGA8_TABLE(EOS).f' 

!*******************************************************************************
      SUBROUTINE SUB_AUTO(PRESS, TEMP, XI) 
!*******************************************************************************      
       IMPLICIT DOUBLE PRECISION (A-Z)     
	   REAL*8 XI(21), TEMP, PRESS     
       RETURN
	   END   	  

	   
!*******************************************************************************	   	
	  SUBROUTINE MOLX(MOL, XIT, XI) 
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: MOLX	  
!*******************************************************************************

!COMPUTES THE MOL WEIGHT, TOTAL OF FRACTIONS AND NORMALIZE IT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL 
      REAL*8  XIT, MCOMP(21)
      MOL=0D0      
      XIT=0D0
      
       MCOMP(1) = 16.04246D0 
       MCOMP(2) = 28.0134D0
       MCOMP(3)= 44.0095D0
       MCOMP(4)= 30.06904D0 
       MCOMP(5)= 44.09562D0
       MCOMP(6)= 58.1222D0
       MCOMP(7)= 58.1222D0
       MCOMP(8)= 72.14878D0
       MCOMP(9)= 72.14878D0
       MCOMP(10)= 86.17536D0
       MCOMP(11)= 100.20194D0
       MCOMP(12)= 114.22852D0
       MCOMP(13)= 128.2551D0
       MCOMP(14)= 142.28168D0
       MCOMP(15)= 2.01588D0
       MCOMP(16)= 31.9988D0
       MCOMP(17)= 28.0101D0
       MCOMP(18)= 18.01528D0
       MCOMP(19)= 34.08088D0
       MCOMP(20)= 4.002602D0
       MCOMP(21)= 39.948D0 
       
      
           
      
!      CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
!     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
!     & COMPR,DENS,D,MOL,A_SUM,GIBBS)    
           
       DO 10 I=1, 21
         XIT = XIT + XI(i)
 10    CONTINUE
       DO 20 I=1,21
         XI(I) =XI(I) * 1d0 / XIT
         MOL = MOL + XI(i)* MCOMP(i)
 20    CONTINUE 
       
      END


!*******************************************************************************
	  SUBROUTINE SOTPX(PRESS, TEMP, ENTRO, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: SOTPX	 
!COMPUTES THE ENTROPY	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL      
	  
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)                    
      END


!*******************************************************************************
	  SUBROUTINE CPOTPX(PRESS, TEMP, CP, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: CPOTPX	
!COMPUTES CP	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D  ,MOL 
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS) 
      END


!*******************************************************************************
	  SUBROUTINE CVOTPX(PRESS, TEMP, CV, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: CVOTPX	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D  ,MOL   
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)              
      END

      
!*******************************************************************************
	  SUBROUTINE WOTPX(PRESS, TEMP, SOUND, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: WOTPX
!COMPUTES SPEED OF SOUND	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D  ,MOL     
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)        
      END


!*******************************************************************************
	  SUBROUTINE CAPVOTPX(PRESS, TEMP, ISENCOEFV, XI)
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: CAPVOTPX
!COMPUTES ISENTROPIC VOLUME COEFFICIENT	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL       
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS) 
      END


!*******************************************************************************
	  SUBROUTINE CAPTOTPX(PRESS, TEMP, ISENCOEFT, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: CAPTOTPX
!COMPUTES ISENTROPIC TEMPERATURE COEFFICIENT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL    
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)
      END


!*******************************************************************************
	  SUBROUTINE RJTOTPX(PRESS, TEMP, JTCOEF, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: RJTOTPX
!COMPUTES JOULE-THOMSON COEFFICIENT	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL      
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)
      END


!*******************************************************************************
	  SUBROUTINE UOTPX(PRESS, TEMP, INTEN, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: UOTPX
!COMPUTES THE INTERNAL ENERGY	(PHYSICAL DIMENSION)	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL       
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)                   
      END


!*******************************************************************************
	  SUBROUTINE HOTPX(PRESS, TEMP, ENTHA, XI)
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: HOTPX
!COMPUTES THE ENTHALPY 	(PHYSICAL DIMENSION)	   
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL       
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)                     
      END


    
!*******************************************************************************
	  SUBROUTINE ZOTPX(PRESS, TEMP, COMPR, XI) 
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: ZOTPX
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: ZOTPX
!COMPUTES THE COMPRESSIBILITY
!USES DL AS COMMON TO CALCULATE COMPRESSIBILITY	  
!INPUT: P,T, XI
!OUTPUT: COMPRESSIBILITY	  
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL  
	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)
      END
      
    
!******************************************************************************* 	  	    
	  SUBROUTINE DOTPX(PRESS, TEMP, DENS, XI) 
!*******************************************************************************
!DEC$ ATTRIBUTES DLLEXPORT, STDCALL, REFERENCE :: DOTPX
!COMPUTES THE DENSITY
!INPUT: P, T, XI
!OUTPUT: DENSITY (DENS) (IN PHYSICAL DIMENSION)	  
      IMPLICIT DOUBLE PRECISION (L-Z)	  
      REAL*8 PRESS, TEMP, XI(21)
	  REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, GIBBS
	  REAL *8 ISENCOEFV, SOUND, COMPR, A_SUM , DENS
      REAL *8 ISENCOEFT,D,MOL        
      
      
 	  CALL AGA8(PRESS, TEMP, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)
      END
      

! Update: 19-apr-2014
! Implement ISO20765 
! Equation of state based on AGA report 8
!==============================================================================
 
        SUBROUTINE AGA8(P, T, XI, INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF,ISENCOEFV,ISENCOEFT,SOUND,
     & COMPR,DENS,D,MOL,A_SUM,GIBBS)
!DEC$ ATTRIBUTES DLLEXPORT:: Gas_Prop
        IMPLICIT DOUBLE PRECISION (A-Z)
       INTEGER i
       REAL*8 XI(21), XJ(21), XJ1(21), T,D,Z, P, BMIX, MOL
       REAL *8 FI0, FI0T, FI0TT, DR, DR_REF, P_REF, T_REF, D_ref
	   REAL *8 INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, ISENCOEFV, SOUND
       REAL *8 FIR, FID, ISENCOEFT, COMPR, A_SUM , DENS, GIBBS
       
     
       DO i=1,21
       XJ1(i) = XI(i)
       ENDDO   
       
       
       XJ(1)=XI(1)
       XJ(2)=XI(2)       
       XJ(3)=XI(3)
       XJ(4)=XI(4)
       XJ(5)=XI(5)
       XJ(6)=XI(18)
       XJ(7)=XI(19)
       XJ(8)=XI(15)
       XJ(9)=XI(17)
       XJ(10)=XI(16)
       XJ(11)=XI(7)
       XJ(12)=XI(6)
       XJ(13)=XI(9)
       XJ(14)=XI(8)       
       XJ(15)=XI(10)       
       XJ(16)=XI(11)
       XJ(17)=XI(12)
       XJ(18)=XI(13)
       XJ(19)=XI(14)
       XJ(20)=XI(20)
       XJ(21)=XI(21)  
                         
       DO i=1,21
       XI(i) = XJ(i)/100d0
       ENDDO 
       
                                       
                       
! =======CALCULATE REFERENCE CONDITIONS======================
  
       	P_REF = 0.101325D0
	    T_REF = 298.15
	    D_REF =0d0
	    DR_REF=0d0
	         
	    CALL DCAGA (XI,MOL, T_REF)
        CALL DZOFPT(P_REF,T_REF,D_REF,Z,BMIX,DR_REF,FIR,FID)
  
!  =========================================================  
    
!CALCULATE HELMOLTZ ENERGY EQUATION PARAMETERS       
        CALL DCAGA (XI,MOL, T)
        CALL DZOFPT(P,T,D,Z,BMIX,DR,FIR,FID)

!CALCULATE HELMOLTZ FREE ENERGY FI_0   
        CALL FI_0(FI0,FI0T,FI0TT,T,T_REF,XI,DR,DR_REF,MOL,D,D_REF) 
!  =========================================================  

!CALCULATE GAS PROPERTIES      
        CALL PROPERTIES(INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, ISENCOEFV, 
     & ISENCOEFT, SOUND, COMPR, DENS, D, T, P, MOL, A_SUM, GIBBS)      
!  =========================================================  

       DO i=1,21
       XI(i) = XJ1(i)
       ENDDO  

      END 


      SUBROUTINE DCAGA (XI,MOL, T)
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER B(58),C(58),K(58),G(58)
      INTEGER Q(58),F(58),S(58),W(58)
      REAL*8 A(58),U(58)
      COMMON /CONSTANTS/ A,B,C,K,U,G,Q,F,S,W
      REAL*8 MW(21),EI(21),KI(21),GI(21),QI(21),FI(21),SI(21),WI(21)
      REAL*8 EIJ(21,21),UIJ(21,21),KIJ(21,21),GIJ(21,21)
      COMMON /PARAMETERS/ MW,EI,KI,GI,QI,FI,SI,WI,EIJ,UIJ,KIJ,GIJ
      REAL*8 K1, CNS(58), CNS1(58), BI(18), BI1(18), BI2(18), BI3(18)
      COMMON /COEF/ K1, CNS, CNS1, BI, BNN, BI1, BI2, BI3
      REAL*8 MWX, RGAS, TCM, DCM, TAU, T
      COMMON /MW/ MWX, RGAS, TCM, DCM
      INTEGER I, J, N
      REAL*8 SUM, XI(21)
      REAL*8 U1, G1, Q1, F1, E1
      REAL*8 XIJ, EIJ0, GIJ0, BNN, MOL
      
      U1=0D0
      G1=0D0
      Q1=0D0 
      F1=0D0 
      E1=0D0
      XIJ=0D0
      EIJ0=0D0
      GIJ0=0D0
      SUM=0D0


      
      DO i=1,21
       DO j=1,21
        IF (EIJ(i,j).EQ.0D0) THEN
         EIJ(i,j)=1D0
        ENDIF
       ENDDO
      ENDDO

      DO i=1,21
       DO j=1,21
        IF (UIJ(i,j).EQ.0D0) THEN
         UIJ(i,j)=1D0
        ENDIF
       ENDDO
      ENDDO
            
      DO i=1,21
       DO j=1,21
        IF (KIJ(i,j).EQ.0D0) THEN
         KIJ(i,j)=1D0
        ENDIF
       ENDDO
      ENDDO
      
      DO i=1,21
       DO j=1,21
        IF (GIJ(i,j).EQ.0D0) THEN
         GIJ(i,j)=1D0
        ENDIF
       ENDDO
      ENDDO      
                           	
      TAU = 1 / T
	
!.....Normalize mole fractions
      SUM = 0
      MWX = 0
      DO 10 I=1, 21
10      SUM = SUM + XI(I)
      DO 20 I=1, 21
20      XI(I) = XI(I)/SUM
!.....Calculate molecular weight
      RGAS = 8.31451D-3
      MWX = 0
      DO 30 I=1, 21
30      MWX = MWX + XI(I)*MW(I) 
      DO 40 N=1, 18
      BI(N) = 0
      BI1(N) =0
      BI2(N) =0
40      BI3(N) =0
      K1 = 0
      U1 = 0
      G1 = 0
      Q1 = 0
      F1 = 0
      E1 = 0
      DO 50 I=1, 21
      K1 = K1 + XI(I)*KI(I)**2.5D0
      U1 = U1 + XI(I)*EI(I)**2.5D0
      G1 = G1 + XI(I)*GI(I)
      Q1 = Q1 + XI(I)*QI(I)
      F1 = F1 + XI(I)*XI(I)*FI(I)
      E1 = E1 + XI(I)*EI(I)
      
50	  CONTINUE
      TCM = 1.261*E1
      DCM = K1**(-1.2D0)
      K1 = K1*K1
      U1 = U1*U1
      DO 60 I=1, 8
      DO 60 J=I+1, 19
      XIJ = XI(I)*XI(J)
      IF (XIJ.NE.0) THEN
      K1 = K1+2.D0*XIJ*(KIJ(I,J)**5.D0-1.D0)*(KI(I)*KI(J))**2.5D0
      U1 = U1+2.D0*XIJ*(UIJ(I,J)**5.D0-1.D0)*(EI(I)*EI(J))**2.5D0
      G1 = G1+XIJ*(GIJ(I,J) - 1.D0)*(GI(I) + GI(J))
      ENDIF
60      CONTINUE
      DO 80 I=1, 21
      DO 80 J=I, 21
      XIJ = XI(I)*XI(J)
      IF (XIJ.NE.0) THEN
      IF (I.NE.J) XIJ = 2.D0*XIJ
      EIJ0 = EIJ(I,J)*DSQRT(EI(I)*EI(J))
      GIJ0 = GIJ(I,J)*(GI(I) + GI(J))/2.D0
      DO 70 N=1, 18
      BNN = (GIJ0 + 1.D0 - G(N))**G(N) 
     &   *(QI(I)*QI(J) + 1.D0 - Q(N))**Q(N)     
     &   *(DSQRT(FI(I)*FI(J)) + 1.D0 - F(N))**F(N)  
     &   *(SI(I)*SI(J) + 1.D0 - S(N))**S(N) 
     &   * (WI(I)*WI(J) + 1.D0 - W(N))**W(N)
      BI(N)=  BI(N)+A(N)*XIJ*EIJ0**U(N)*(KI(I)*
     &   KI(J))**1.5D0*BNN*TAU**U(N)
        BI1(N)=BI1(N)+U(N)*A(N)*XIJ*EIJ0**U(N)*
     &   (KI(I)*KI(J))**1.5D0*BNN*TAU**U(N)
		BI2(N) = BI2(N)+(U(N)**2-U(N))*A(N)*XIJ*EIJ0**U(N)*
     &   (KI(I)*KI(J))**1.5D0*BNN*TAU**U(N)
		BI3(N) = BI3(N)+(1-U(N))*A(N)*XIJ*EIJ0**U(N)*(KI(I)
     &   *KI(J))**1.5D0*BNN*TAU**U(N)       

70      CONTINUE
      ENDIF
80      CONTINUE
     
      K1 = K1**0.2D0	
      U1 = U1**0.2D0
      DO 90 N=13, 58
90      CNS(N) = (G1 + 1.D0 - G(N))**G(N)  
     &   * (Q1**2 + 1.D0 - Q(N))**Q(N) 
     &   * (F1 + 1.D0 - F(N))**F(N)
     &   * A(N)*U1**U(N)

      MOL = MWX    	  
      END
!======================================================================


	   SUBROUTINE PZOFDT(D, T, P, Z,BMIX,DR,FIR,FID)
	   IMPLICIT DOUBLE PRECISION (A-Z)
	   INTEGER B(58),C(58),K(58),G(58)
	   INTEGER Q(58),F(58),S(58),W(58)

	   REAL*8 A(58),U(58)

	   COMMON /CONSTANTS/ A,B,C,K,U,G,Q,F,S,W
	   REAL*8 K1, CNS(58), CNS1(58), BI(18), BI1(18),BI2(18),BI3(18)
	   REAL*8 BNN,TAU
	   REAL*8 FIR,FID, FI_1, FIR_T, FIR_TT, FI_2
	   COMMON /COEF/ K1, CNS, CNS1, BI, BNN, BI1, BI2, BI3
	   REAL*8 MWX, RGAS, TCM, DCM, DR
	   REAL*8 DERI_HELM(9)
	   COMMON /MW/ MWX, RGAS, TCM, DCM
       COMMON /DERIVATIVES / DERI_HELM
    	
	   INTEGER N
	   REAL*8 D, T, P, Z, BMIX
	   fir =0d0
       TAU = 1 / T 
	
	   DR = D*K1**3
	   BMIX = 0
	   DO 10 N=1, 18
10	   BMIX = BMIX + BI(N)

       Z = 1.D0 + BMIX*D! *TAU**U(N)
	  
       DO 20 N=13, 18
20	     Z = Z - DR*CNS(N)*TAU**U(N)
       DO 30 N=13, 58
	 
30	     Z = Z + CNS(N)*TAU**U(N)*(B(N) - 
     &   C(N)*K(N)*DR**K(N))*DR**B(N)
     &   *DEXP(-C(N)*DR**K(N))
       P = D*RGAS*Z / TAU

! !THIS CALCULATE FIR (helmholtz residual energy)	
       DO 50 N=13, 18
       FIR = FIR - DR*CNS(N)*TAU**U(N)
50      continue

       DO 60 N=13, 58
60       FIR = FIR + CNS(N)*TAU**U(N)*DR**B(N)  
     &   *DEXP(-C(N)*DR**K(N))    
       FIR = FIR +  BMIX * D

! !THIS CALCULATE FID
       FID = 0D0
       FID = BMIX * D 
       DO 70 N=13, 18
70       FID = FID - DR*CNS(N)*TAU**U(N)

       DO 80 N=13, 58
80	     FID = FID + CNS(N)*TAU**U(N)*(B(N) - 
     &   C(N)*K(N)*DR**K(N))*DR**B(N) 
     &   *DEXP(-C(N)*DR**K(N))
     
       FID = FID
      
! !THIS CALCULATE FI_1
 
       FI_1 = 2 * BMIX * D + 1D0
              
       DO 90 N=13, 18
90	   FI_1 = FI_1 - 2 * DR*CNS(N)*TAU**U(N)

       DO 100 N=13, 58
100	    FI_1 = FI_1 + CNS(N)*TAU**U(N)*
     &   ( B(N)-(1+K(N))*C(N)*K(N)*DR**K(N)+ 
     &   (B(N) - C(N)*K(N)*DR**K(N))**2  ) *DR**B(N)
     &   *DEXP(-C(N)*DR**K(N))
        
!THIS CALCULATE FIR_T
        
       FIR_T = 0D0
        
       DO 110 N=1, 18
110	   FIR_T = FIR_T + D * BI1(N)


       DO 115  N=13, 18
115	    FIR_T  = FIR_T  - DR * U(N)*CNS(N)*TAU**U(N)     
       

       DO 120 N=13, 58
120 	 FIR_T = FIR_T + U(N) * CNS(N)*(TAU**U(N))*(DR**B(N))
     &   *DEXP(-C(N)*DR**K(N))

!THIS CALCULATE FIR_TT
       FIR_TT = 0D0
        
       DO 130 N=1, 18
130	   FIR_TT = FIR_TT + D * BI2(N)

 
   	   DO 140  N=13, 18
140	   FIR_TT  = FIR_TT  - DR * (U(N)**2-U(N))*CNS(N)*TAU**U(N)     
       

       DO 150 N=13, 58
150 	 FIR_TT = FIR_TT + (U(N)**2-U(N))* CNS(N)*(TAU**U(N))*(DR**B(N))
     &   *DEXP(-C(N)*DR**K(N))
        
!THIS CALCULATE FI2
       FI_2 = 1D0       
       DO 160 N=1, 18
160	     FI_2 = FI_2  + D * BI3(N)
  
       DO 170  N=13, 18
170	    FI_2  = FI_2  - DR * (1-U(N))*CNS(N)*TAU**U(N)     
       
       DO 180 N=13, 58
180 	 FI_2 = FI_2 + (1-U(N))* CNS(N)*(TAU**U(N))*(DR**B(N))
     &   *(B(N) - C(N)*K(N)*DR**K(N))*DEXP(-C(N)*DR**K(N))
       
       DERI_HELM(1) = FIR_T
       DERI_HELM(2) = FIR_TT
       DERI_HELM(3) = FID
       DERI_HELM(4) = FI_1
       DERI_HELM(5) = FI_2
       DERI_HELM(6) = FIR

	END

!======================================================================
	   SUBROUTINE DZOFPT(P, T, D, Z, BMIX, DR, FIR, FID)
	   IMPLICIT DOUBLE PRECISION (A-Z)
	   REAL*8 P, T, D, Z, BMIX, DR
	   REAL*8 X1, X2, X3, F, F1, F2, F3, TOL
	   REAL*8 FIR, FID
	   INTEGER I
	
	   TOL = 0.5D-9
	   X1 = 0.000001D0
	   X2 = 40.D0
	   D = 0
	   CALL PZOFDT(X1, T, F1, Z, BMIX, DR, FIR, FID)
	   CALL PZOFDT(X2, T, F2, Z, BMIX, DR, FIR, FID)
	   F1 = F1 - P
	   F2 = F2 - P
	   IF (F1*F2.GE.0) RETURN
!----------------------------------------------------------------------
!BEGIN ITERATING
!----------------------------------------------------------------------
	  DO 60 I = 1, 50
!...Use False Position to get point 3.
	  X3 = X1 - F1*(X2 - X1)/(F2 - F1)
	  CALL PZOFDT(X3, T, F3, Z, BMIX, DR, FIR, FID)
	  F3 = F3 - P
!...Use points 1, 2, and 3 to estimate the root using Chamber's
!...method (quadratic solution).
	  D = X1*F2*F3/((F1 - F2)*(F1 - F3)) 
     &   + X2*F1*F3/((F2 - F1)*(F2 - F3)) 
     &   + X3*F1*F2/((F3 - F1)*(F3 - F2))
	  IF ((D - X1)*(D - X2).GE.0) D = (X1 + X2)/2.D0
	  CALL PZOFDT(D, T, F, Z, BMIX, DR, FIR, FID)
	  F = F - P
	  IF (DABS(F).LE.TOL) RETURN
!...Discard quadratic solution if false position root is closer.
	  IF (DABS(F3).LT.DABS(F) .AND. F*F3.GT.0) THEN
	  IF (F3*F1.GT.0) THEN
	  X1 = X3
	  F1 = F3
	  ELSE
	  X2 = X3
	  F2 = F3
	  ENDIF
	  ELSE
!...Swap in new value from quadratic solution
	  IF (F*F3.LT.0) THEN
	  X1 = D
	  F1 = F
	  X2 = X3
	  F2 = F3
	  ELSEIF (F3*F1.GT.0) THEN
	  X1 = D
	  F1 = F
	  ELSE
	  X2 = D
	  F2 = F
	  ENDIF
	  ENDIF
60	  CONTINUE
	  D = 0
	  END

      SUBROUTINE FI_0(FI0,FI0T,FI0TT,T,T_REF,XI,DR,DR_REF,MOL,D,D_REF)
      IMPLICIT DOUBLE PRECISION (A-Z)      
      REAL*8  FI0, FI0T, FI0TT, TAU, T, DR, DR_REF, TAU_REF
      REAL*8  MOL, D, T_REF, D_REF
      REAL*8 A01(21),A02(21),B0(21), C0(21), D0(21), E0(21)
      REAL*8 XI(21), DERI_HELM(9)  
      REAL*8 F0(21), G0(21), H0(21), I0(21), J0(21) 
      COMMON /PARAMETERSFREE/ A01,A02,B0, C0, D0, E0, 
     &   F0, G0, H0, I0, J0
      COMMON /DERIVATIVES / DERI_HELM       
      INTEGER I

      TAU = 1 / T
      TAU_REF = 1/ T_REF	     
	     
!THIS CALCULATE FI0  (helmholtz free energy)   	        
      FI0 = 0D0

      DO 200 i=1,21          
       IF (XI(i).eq.0d0) THEN
       GOTO 195
      ENDIF

      FI0 = FI0 + XI(i)*(A01(i) + A02(i) * TAU + B0(i) * log(tau) 
     &   +C0(i) * Log(Sinh(D0(i)*TAU)) - E0(i) * LOG(Cosh(F0(i)*Tau)) 
     &   +G0(i) * Log(Sinh(H0(i)*Tau)) - I0(i) * Log(Cosh(J0(i)*Tau)) 
     &   +Log(XI(i)) ) 	      
	  

195      CONTINUE
200      CONTINUE

      FI0 = FI0 + Log(DR/DR_REF)  + Log(TAU_REF/ TAU)
      
!THIS CALCULATE FI0T       
      FI0T = 0D0
                  
      DO 210 i=1,21

      FI0T = FI0T + XI(i)*( A02(i) + (B0(i)-1) / tau 
     &   +  C0(i) * D0(i) * ( cosh(D0(i)*TAU) / Sinh(D0(i)*TAU)  ) 
     &   -  E0(i) * F0(i) * ( sinh(F0(i)*Tau) / cosh(F0(i)*Tau)  ) 
     &   +  G0(i) * H0(i) * ( cosh(H0(i)*TAU) / Sinh(H0(i)*TAU)  )
     &   -  I0(i) * J0(i) * ( sinh(J0(i)*Tau) / cosh(J0(i)*Tau) )) 
210      CONTINUE
         
      FI0T =   FI0T * tau
          
!THIS CALCULATE FI0TT  
      FI0TT = 0D0
              
      DO 220 i=1,21
	      FI0TT = FI0TT + XI(i)*( -(B0(i)-1) / tau**2 
     &   - C0(i) *( D0(i) / Sinh(D0(i)*TAU)  )**2    
     &   - E0(i) *( F0(i) / Cosh(F0(i)*TAU)  )**2    
     &   - G0(i) *( H0(i) / Sinh(H0(i)*TAU)  )**2    
     &   - I0(i) *( J0(i) / Cosh(J0(i)*TAU)  )**2  )  
220      CONTINUE
         
      FI0TT =   FI0TT * tau**2

      DERI_HELM(7) = FI0
      DERI_HELM(8) = FI0T
      DERI_HELM(9) = FI0TT
    
      END
	
!=======================================================================	
       SUBROUTINE PROPERTIES(INTEN, ENTHA, ENTRO, CV, CP, 
     & JTCOEF, ISENCOEFV, ISENCOEFT, SOUND, COMPR, DENS, 
     & D, T, P, MOL, A_SUM, GIBBS)   
       IMPLICIT DOUBLE PRECISION (A-Z)       
       REAL*8 ISENCOEFV,SOUND,COMPR,T,MOL,D
       REAL*8 A_SUM,ISENCOEFT, GIBBS
       REAL*8 DERI_HELM(9), DENS,INTEN, ENTHA, ENTRO, CV, CP, JTCOEF, P
      
       COMMON /DERIVATIVES / DERI_HELM  

!      for reference:
!      DERI_HELM(1) = FIR_T
!      DERI_HELM(2) = FIR_TT
!      DERI_HELM(3) = FID
!      DERI_HELM(4) = FI_1
!      DERI_HELM(5) = FI_2
!	   DERI_HELM(6) = FIR
!      DERI_HELM(7) = FI0
!      DERI_HELM(8) = FI0T
!      DERI_HELM(9) = FI0TT

      RGAS = 8.31451D-3

      COMPR = 1 + DERI_HELM(3)
 
      GIBBS = (COMPR+DERI_HELM(6)+DERI_HELM(7))* ((RGAS *1000)* T / MOL)
      
      INTEN = (DERI_HELM(8) + DERI_HELM(1))  
     &   * ((RGAS *1000) * T / MOL)

      ENTHA = INTEN  + COMPR * (RGAS *1000) * T / MOL
 
      ENTRO=  INTEN / T -  (DERI_HELM(7) + DERI_HELM(6)) 
     &   * ((RGAS *1000) / MOL)

      CV = (RGAS *1000) *  
     &   (-DERI_HELM(9)- DERI_HELM(2))  / MOL

      CP =    CV + (DERI_HELM(5)**2/ DERI_HELM(4) )
     &   *  (RGAS *1000)  / MOL
		
      JTCOEF =  (1000 / (CP*D*MOL)) 
     &   * ((DERI_HELM(5)/ DERI_HELM(4))- 1)
		
      ISENCOEFV = (DERI_HELM(4)/COMPR) * CP/CV

      SOUND =  SQRT(DERI_HELM(4) * 
     &   ((RGAS *1000) * T / MOL) * CP *1000/ (CV))

      A_SUM = (DERI_HELM(7) + DERI_HELM(6))* 
     &   ((RGAS *1000) * T / MOL)
     
      ISENCOEFT=T / (T-1D3*P/(D* MOL*CP)-P*JTCOEF)
      
      DENS = D * MOL
      END





	  
