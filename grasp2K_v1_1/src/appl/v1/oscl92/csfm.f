************************************************************************
*                                                                      *
      SUBROUTINE CSFM (ASFA,ASFB,LEV1,LEV2)
*                                                                      *
*   This routine calculates  the CSF Coulomb, Babuskin, and magnetic   *
*   matrix elements for  a transition  between  levels  separated by   *
*   energy OMEGA.                                                      *
*                                                                      *
*   Call(s) to: [OSCL92]: SPME.                                        *
*               [LIB92]: ITJPO.                                        *
*                                         Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNTJJA,JJA(1))
      POINTER (PNTJJB,JJB(1))
      POINTER (PNTHB1,HB1(1))
      POINTER (PNTHC1,HB2(1))
      POINTER (PNTHB2,HC1(1))
      POINTER (PNTHC2,HC2(1))
      POINTER (PNTHM1,HM1(1))
      POINTER (PNTHM2,HM2(1))
      POINTER (PISLDR,ISLDR(1))
      POINTER (PXSLDR,XSLDR(1))
      POINTER (PNTLAB,LAB(1))
      POINTER (PNNPTR,NPTR(1))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /OSC1/PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :            PNTHC1,PNTHC2,PNTHM1,PNTHM2,NSDIM
     :      /OSC2/LK,KK
     :      /OSC3/PXSLDR,PISLDR,NTDIM
     :      /OSC5/NINT,PNTLAB,PNNPTR,NINTEG
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB,
     :           NCA = 65536)
*
      EPS = 1.0D-10
*
      ASFA = 0.0D 00
      ASFB = 0.0D 00
      NCOUNT = 0
*
      DO 4 M = 1,NINTEG
         JA = LAB(M)
         JB = MOD (JA,KEY)
         JA = JA/KEY
         NCLR = NPTR(M)+1
         NCUP = NPTR(M+1)
         IF (NCLR .GT. NCUP) GOTO 4
         IF (MOD(NKL(JA)+NKL(JB)+KK+LK,2) .NE. 0) GOTO 4
         IF (JA .LT. JB) THEN
            NCOUNT = NCOUNT+1
            JJA(NCOUNT) = JA
            JJB(NCOUNT) = JB
            CALL SPME (JA,JB,HCOUL1,HBAB1,HMAG1)
            HC1(NCOUNT) = HCOUL1
            HB1(NCOUNT) = HBAB1
            HM1(NCOUNT) = HMAG1
            CALL SPME (JB,JA,HCOUL2,HBAB2,HMAG2)
            HC2(NCOUNT) = HCOUL2
            HB2(NCOUNT) = HBAB2
            HM2(NCOUNT) = HMAG2
         ELSEIF (JA .EQ. JB) THEN
            CALL SPME (JA,JB,HCOUL1,HBAB1,HMAG1)
            HCOUL2 = HCOUL1
            HBAB2 = HBAB1
            HMAG2 = HMAG1
         ELSE
            DO 1 J = 1,NCOUNT
               IF ((JA .EQ. JJB(J)) .AND. (JB .EQ. JJA(J))) THEN
                  HCOUL1 = HC2(J)
                  HBAB1 = HB2(J)
                  HMAG1 = HM2(J)
                  HCOUL2 = HC1(J)
                  HBAB2 = HB1(J)
                  HMAG2 = HM1(J)
                  GOTO 2
               ENDIF
    1       CONTINUE
            CALL SPME (JA,JB,HCOUL1,HBAB1,HMAG1)
            CALL SPME (JB,JA,HCOUL2,HBAB2,HMAG2)
         ENDIF
    2    CONTINUE
*
         DO 3 J = NCLR,NCUP
*
            IA = ISLDR(J)
CFF         .. note that if the expansion is too long some
*           packed indices will be treated as negative numbers
*           Could this extension be architecture dependent?
	    IF (IA .GT. 0) THEN
              IB = MOD (IA,NCA)
              IA = IA/NCA
	    ELSE
	      IB = NCA + MOD(IA, NCA)
	      IA = NCA + IA/NCA - 1
	    END IF
*
            COUVX = EVEC(IA+(LEV1-1)*NCF)*EVEC(IB+(LEV2-1)*NCF)
            COEFF = XSLDR(J)
            IF (ABS(COUVX) .GT. EPS) THEN
               IF (LDBPR(18)) WRITE (99,300) NP(JA),NH(JA),
     :                                       NP(JB),NH(JB),
     :                                       IA,IB,COEFF
               IF (KK .EQ. 0) THEN
                  ASFA = ASFA+HCOUL1*COEFF*COUVX
                  ASFB = ASFB+HBAB1*COEFF*COUVX
               ELSE
                  ASFA = ASFA+HMAG1*COEFF*COUVX
               ENDIF
            ENDIF
            COUVX = EVEC(IB+(LEV1-1)*NCF)*EVEC(IA+(LEV2-1)*NCF)
            IF ((ABS (COUVX) .GT. EPS) .AND. (IA .NE. IB)) THEN
               FACT1 = ITJPO(IA)*(NKJ(JB)+1)
               FACT2 = ITJPO(IB)*(NKJ(JA)+1)
               FACT = SQRT(FACT1/FACT2)
               IDL = (ITJPO(IB)-ITJPO(IA)+NKJ(JB)-NKJ(JA))/2
               IF (MOD(IDL,2) .NE. 0) FACT = -FACT
               COEFF = FACT*COEFF
               IF (LDBPR(18)) WRITE (99,300) NP(JB),NH(JB),
     :                                       NP(JA),NH(JA),
     :                                       IB,IA,COEFF
               IF (KK .EQ. 0) THEN
                  ASFA = ASFA+HCOUL2*COEFF*COUVX
                  ASFB = ASFB+HBAB2*COEFF*COUVX
               ELSE
                  ASFA = ASFA+HMAG2*COEFF*COUVX
               ENDIF
            ENDIF
    3    CONTINUE
    4 CONTINUE
*
      RETURN
*
  300 FORMAT ('MCT coefficient',2(1X,1I2,1A2),2I8,1P,D19.8)
*
      END

