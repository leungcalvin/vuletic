************************************************************************
*                                                                      *
      SUBROUTINE BREIT (JA,JB,JA1,JB1,JA2,JB2)
*                                                                      *
*   Computes  the  coefficients  appearing in EQS. 5, 8, 9 AND 10 OF   *
*   I P Grant and B J McKenzie,  J Phys B 13 (1980) 2671--2681.  The   *
*   coefficients for each choice of  orbitals JA1, JB1, JA2, and JB2   *
*   depend on two further parameters  NU and K;  there are IMU inte-   *
*   grals for each such choice, where:                                 *
*                                                                      *
*                  IMU = 4          TYPE = 1                           *
*                        8                 2                           *
*                        1                 3                           *
*                        1                 4                           *
*                        3                 5                           *
*                        4                 6                           *
*                                                                      *
*   See the paper cited above for details.                             *
*                                                                      *
*   Call(s) to: [LIB92]: GENSUM, ITRIG, KNJ, LTAB, MODJ23, MUMDAD,     *
*                        NJGRAF, OCON, SETJ                            *
*               [RCI92]: CXK, SNRC, TALK.                              *
*                                                                      *
*                                           Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
*
      PARAMETER (
     :   MANGM=60,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3),
     :   MTRIAD=12,M2TRD=2*MTRIAD,M4TRD=4*MTRIAD,
     :   M6J=20,MSUM=10)
*
      LOGICAL FREE,DSUMVR,ESUMVR,FAILD,FAILE
*
      DIMENSION COND(12,20),CONE(12,20),S(12)
      DIMENSION IS(4),KAPS(4),KS(4),NQS(4),ILS(4),LLS(4),IT1(4),IROWS(4)
      DIMENSION JS(4)
*
      DIMENSION JD6(M3MNGM),JD7(M3MNGM),JD8(M3MNGM),
     : JD9(MANGMP),KDW(6,M6J),LDDEL(M6J,2),DSUMVR(MANGM)
      DIMENSION JD6P(MANGMP),JD7P(MANGMP),JD8P(MANGMP),JD9P(MANGMP),
     : JDWORD(6,M6J),
     : NDBJ(MSUM),NDB6J(MSUM),KD6CP(MSUM),KD7CP(MSUM),KD8CP(MSUM),
     : KD9CP(MSUM),JDSUM6(MTRIAD),JDSUM4(MTRIAD,M6J),JDSUM5(MTRIAD,M6J),
     : INVD6J(M6J)
      DIMENSION JE6(M3MNGM),JE7(M3MNGM),JE8(M3MNGM),
     : JE9(MANGMP),KEW(6,M6J),LEDEL(M6J,2),ESUMVR(MANGM)
      DIMENSION JE6P(MANGMP),JE7P(MANGMP),JE8P(MANGMP),JE9P(MANGMP),
     : JEWORD(6,M6J),
     : NEBJ(MSUM),NEB6J(MSUM),KE6CP(MSUM),KE7CP(MSUM),KE8CP(MSUM),
     : KE9CP(MSUM),JESUM6(MTRIAD),JESUM4(MTRIAD,M6J),JESUM5(MTRIAD,M6J),
     : INVE6J(M6J)
*
      COMMON/COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /L1/JBQ1(3,NNNW),JBQ2(3,NNNW),JTQ1(3),JTQ2(3)
     :      /L2/J2S(MTRIAD,3),J3S(MTRIAD,3)
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      PARAMETER (EPS = 1.0D-10)
      PARAMETER (NUMAX = 20)
*
*   1.0  Initialize pointers and flags and set any
*        tables required.
*
*        In this segment, the array IS points to the
*        full list of orbitals, the array JS to the
*        array JLIST of peel orbital pointerS.
*
*   1.1  Initialization
*
      JS(1) = JA1
      JS(2) = JB1
      JS(3) = JA2
      JS(4) = JB2
      DO 1 I = 1,4
        IS(I) = JLIST(JS(I))
        KAPS(I) = 2*NAK(IS(I))
        KS(I) = ABS (KAPS(I))
    1 CONTINUE
      IA1 = IS(1)
      IB1 = IS(2)
      IA2 = IS(3)
      IB2 = IS(4)
      NQS(1) = NQ1(IA1)
      NQS(2) = NQ1(IB1)
      NQS(3) = NQ2(IA2)
      NQS(4) = NQ2(IB2)
*
      KJ23 = 0
      ISNJ = 0
*
      FAILD = .FALSE.
      FAILE = .FALSE.
      NBRJ = 3*NPEEL + 7
      DO 2 I = 1,(NBRJ-1)
        FREE(I) = .FALSE.
    2 CONTINUE
      FREE(NBRJ) = .TRUE.
*
*   2.0  Set quantum numbers of spectator shells.
*
      DO 4 J = 1,NW
        DO 3 K = 1,3
          JBQ1(K,J) = 0
          JBQ2(K,J) = 0
    3   CONTINUE
    4 CONTINUE
*
      DO 8 JJ = 1,NPEEL
        J = JLIST(JJ)
        IF ((J .NE. IA1) .AND. (J .NE. IB1)) THEN
          DO 5 K = 1,3
            JBQ1(K,J) = JJQ1(K,J)
    5     CONTINUE
        ENDIF
        IF ((J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 6 K = 1,3
            JBQ2(K,J) = JJQ2(K,J)
    6     CONTINUE
        ENDIF
*
*   2.1  Examine spectator shells for orthogonality
*
        IF ((J .NE. IA1) .AND. (J .NE. IB1) .AND.
     :      (J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 7 K = 1,3
            IF (JBQ1(K,J) .NE. JBQ2(K,J) ) GOTO  98
    7     CONTINUE
        ENDIF
    8 CONTINUE
*
*   3.0  Start main calculation
*        Begin with common factors
*
      CONST = OCON (IA1,IB1,IA2,IB2)
      IF (IBUG2 .NE. 0) WRITE (99,307) CONST
*
*   3.1  Set range of tensor index NU
*
      IF (IBUG2 .NE. 0) WRITE (99,302) IA1,IB1,IA2,IB2
      CALL SNRC (IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
      IF (IBUG2 .NE. 0) WRITE (99,303) ND1,ND2,NE1,NE2,IBRD,IBRE
      IF ((IBRD .LT. 0) .AND. (IBRE .LT. 0)) RETURN
      IF ((ND2 .GT. NUMAX) .OR. (NE2 .GT. NUMAX)) THEN
         KK = MAX (ND2,NE2)
         WRITE (*,301) KK
         STOP
      ENDIF
      IF (IBRD .GE. 0) THEN
        DO 10 N = 1,ND2
          DO 9 MU = 1,12
            COND(MU,N) = 0.0D 00
    9     CONTINUE
   10   CONTINUE
      ENDIF
      IF (IBRE .GE. 0) THEN
        DO 12 N = 1,NE2
          DO 11 MU = 1,12
            CONE(MU,N) = 0.0D 00
   11     CONTINUE
   12   CONTINUE
      ENDIF
*
*   3.2  Set parameters of summation over parent
*        (barred) terms in Eq. 2 (loc cit). The array
*        IROWS is formed to point to the list of
*        allowed parents of active shells in the
*        array NTAB
*
      CALL LTAB (IS,NQS,KS,IROWS)
*
      DO 13 I = 1,4
        II = IROWS(I)
        LLS(I) = ITAB(II)
        ILS(I) = JTAB(II)
   13 CONTINUE
*
*   4.0  Sum over all parent terms permitted by
*        angular momentum and seniority selection rules
*
      LLS1 = LLS(1)
      IF (LLS1 .NE. 1) FREE(JA1) = .TRUE.
      LLS2 = LLS(2)
      LLS3 = LLS(3)
      LLS4 = LLS(4)
*
      LS2 = ILS(2)
      DO 29 LB1 = 1,LLS2
        LS2 = LS2+3
        IT1(2) = NTAB(LS2)
        IT12 = IT1(2)
        IT2 = KS(2)
        IT3 = JJQ1(3,IB1)
        IF (ITRIG (IT12,IT2,IT3) .EQ. 0) GOTO 29
        IF (ABS (NTAB(LS2-2)-JJQ1(1,IB1) ) .NE. 1) GOTO 29
*
        LS1 = ILS(1)
        DO 28 LA1 = 1,LLS1
          LS1 = LS1+3
          IT1(1) = NTAB(LS1)
          IT11 = IT1(1)
          IT2 = KS(1)
          IF (IA1 .EQ. IB1) THEN
*
*   Treat IA1 .EQ. IB1 as a special case
*
            IT3 = IT1(2)
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 28
            IF (ABS (NTAB(LS1-2)-NTAB(LS2-2)) .NE. 1) GOTO 28
            IF (LLS2 .NE. 1) FREE(NBRJ-8) = .TRUE.
          ELSE
            IT3 = JJQ1(3,IA1)
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 28
            IF (ABS (NTAB(LS1-2)-JJQ1(1,IA1)) .NE. 1) GOTO 28
            IF (LLS2 .NE. 1) FREE(JB1) = .TRUE.
          ENDIF
*
          LS4 = ILS(4)
          DO 27 LB2 = 1,LLS4
            LS4 = LS4+3
            IT1(4) = NTAB(LS4)
            IT14 = IT1(4)
            IT2 = KS(4)
            IT3 = JJQ2(3,IB2)
            IF (ITRIG(IT14,IT2,IT3) .EQ. 0) GOTO 27
            IF (ABS (NTAB(LS4-2)-JJQ2(1,IB2)) .NE. 1) GOTO 27
*
            LS3 = ILS(3)
            DO 26 LA2 = 1,LLS3
              LS3 = LS3+3
              IT1(3) = NTAB(LS3)
              IT13 = IT1(3)
              IT2 = KS(3)
              IF (IA2 .EQ. IB2) THEN
*
*   TREAT IA2 .EQ. IB2 as a special case
*
                IT3 = IT1(4)
                IF (LLS4 .NE. 1) FREE(NBRJ-6) = .TRUE.
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 26
                IF (ABS (NTAB(LS3-2)-NTAB(LS4-2)) .NE. 1) GOTO 26
              ELSE
                IT3 = JJQ2(3,IA2)
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 26
                IF (ABS (NTAB(LS3-2)-JJQ2(1,IA2)) .NE. 1) GOTO 26
              ENDIF
*
*   At this point the current parent has been completely defined,
*   and its quantum numbers can now be set.  The JTQ arrays must
*   be set if IA1 .EQ. IB1 or IA2 .EQ. IB2. The matrix element should be
*   diagonal in barred quantum numbers.
*
              DO 15 K = 1,3
                JBQ1(K,IA1) = NTAB(LS1+K-3)
                JBQ2(K,IA2) = NTAB(LS3+K-3)
                JTQ1(K) = 0
                IF (IB1 .EQ. IA1) THEN
                  JTQ1(K) = NTAB(LS2+K-3)
                ELSE
                  JBQ1(K,IB1) = NTAB(LS2+K-3)
                ENDIF
                JTQ2(K) = 0
                IF (IB2 .EQ. IA2) THEN
                  JTQ2(K) = NTAB(LS4+K-3)
                ELSE
                  JBQ2(K,IB2) = NTAB(LS4+K-3)
                ENDIF
                DO 14 KK = 1,4
                  IF (JBQ1(K,IS(KK)) .NE. JBQ2(K,IS(KK))) GOTO 26
   14           CONTINUE
   15         CONTINUE
*
*   4.1 Evaluate product of 4 CFPs
*
              CALL MUMDAD (IS,KAPS,PROD)
              IF (ABS (PROD) .LT. EPS) GOTO 26
*
*    4.2  Set arrays for defining the recoupling
*         coefficient
*
              CALL SETJ (IS,JS,KS,NPEEL,KJ23)
*
              IF (ISNJ .EQ. 0) THEN
*
********************* N J G R A F   V E R S I O N **********************
*
*     Set up the arrays and variables for the direct case.
*
                IF (IBRD.GE.0) THEN
                  CALL NJGRAF (RECUP,-1,FAILD)
                  ISNJ = 1
                  IF (.NOT. FAILD) THEN
                  CALL KNJ(JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,
     :                     KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                     JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,
     :                     NDB6J,KD6CP,KD7CP,KD8CP,KD9CP,
     :                     JDSUM4,JDSUM5,JDSUM6,INVD6J)
                  ENDIF
                ENDIF
*
*   Set up the arrays and variables for the exchange case.
*
                IF (IBRE .GE. 0) THEN
                  CALL MODJ23
                  CALL NJGRAF (RECUP,-1,FAILE)
                  ISNJ = 2
                  IF (.NOT. FAILE) THEN
                  CALL KNJ(JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,JE9,
     :                     KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                     JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,
     :                     NEB6J,KE6CP,KE7CP,KE8CP,KE9CP,
     :                     JESUM4,JESUM5,JESUM6,INVE6J)
                  ENDIF
                ENDIF
*
              IF (FAILD .AND. FAILE) GOTO 30
*
              ENDIF
*
*   4.3.1 Summation for direct terms
*
              IF ((IBRD .GE. 0). AND. (.NOT. FAILD)) THEN
*
                IMUD = 4
                IF (IBRD .GT. 1) IMUD = 1
                NCODE = 0
                DO 20 N = 1,ND2
                  NU = ND1+2*(N-1)
                  NUD = NU+NU+1
*
                  IF (NU .NE. 0) THEN
*
                    IF ((ITRIG (KS(1),KS(3),NUD) .NE. 0) .AND.
     :                  (ITRIG (KS(2),KS(4),NUD) .NE. 0)) THEN
*
                    K = NU
                    J1(MJA) = NUD
                    CALL GENSUM (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,
     :                JD9,KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :                KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,
     :                INVD6J,X)
                    IF (IBUG2 .NE. 0) WRITE (99,304) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK(S,IS,KAPS,NU,K,IBRD,1)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,IMUD)
                      DO 16 MU = 1,IMUD
                        COND(MU,N) = COND(MU,N)+X*S(MU)
   16                 CONTINUE
                    ENDIF
*
                  ENDIF
*
*   K = NU-1
*
                  IF (IBRD .GT. 1) GOTO 20
*
                  K = NU-1
*
                  IF (NCODE .EQ. N) THEN
                    X=XCODE
                  ELSE
                    ITKMO = NUD-2
                    IF (ITRIG (KS(1),KS(3),ITKMO) .EQ. 0) GOTO 18
                    IF (ITRIG (KS(2),KS(4),ITKMO) .EQ. 0) GOTO 18
                    J1(MJA) = ITKMO
                    CALL GENSUM(JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,
     :                JD9,KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :                KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,
     :                INVD6J,X)
                  ENDIF
*
                  IF (IBUG2 .NE. 0) WRITE (99,304) NU,K,X
*
                  IF (ABS (X) .GE. EPS) THEN
                    X = X*PROD
                    CALL CXK (S,IS,KAPS,NU,K,IBRD,1)
                    IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                (S(III),III = 1,4)
                    DO 17 MU = 1,4
                      COND (MU,N) = COND(MU,N)+X*S(MU)
   17               CONTINUE
                  ENDIF
*
                ENDIF
*
*   K = NU+1
*
   18           IF ((IBRD .GT. 1) .OR. (N .EQ. ND2)) GOTO 20
*
                NCODE = N+1
                XCODE = 0.0D 00
                ITKMO = NUD+2
*
                IF ((ITRIG(KS(1),KS(3),ITKMO) .NE. 0) .AND.
     :              (ITRIG(KS(2),KS(4),ITKMO) .NE. 0)) THEN
                  K = NU+1
                  J1(MJA) = ITKMO
                  CALL GENSUM (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,
     :              JD9,KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :              JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :              KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,
     :              INVD6J,X)
                    XCODE = X
                    IF (IBUG2 .NE. 0) WRITE (99,304) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK (S,IS,KAPS,NU,K,IBRD,1)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,12)
                      DO 19 MU = 1,12
                        COND(MU,N) = COND(MU,N)+X*S(MU)
   19                 CONTINUE
                    ENDIF
*
                  ENDIF
*
   20           CONTINUE
*
              ENDIF
*
*   4.3.2 Summation for exchange terms
*
              IF ((IBRE .GE. 0) .AND. (.NOT. FAILE)) THEN
*
                NCODE = 0
*
                DO 25 N = 1,NE2
                  IMUE = 4
                  IF (IBRE .EQ. 2) IMUE = 1
                  IF (IBRE .EQ. 4) IMUE = 3
                  NU = NE1+2*(N-1)
                  NUD = NU+NU+1
*
                  IF (NU .NE. 0) THEN
*
                    IF ((ITRIG(KS(1),KS(4),NUD) .NE. 0) .AND.
     :                  (ITRIG(KS(2),KS(3),NUD) .NE. 0)) THEN
                      K = NU
                      J1(MJA) = NUD
                      CALL GENSUM(JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,
     :                  JE9,KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                  JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     :                  KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,
     :                  INVE6J,X)
                      IF (IBUG2 .NE. 0) WRITE (99,306) NU,K,X
*
                      IF (ABS (X) .GE. EPS) THEN
                        X = X*PROD
                        CALL CXK (S,IS,KAPS,NU,K,IBRE,2)
                        IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                    (S(III),III = 1,IMUE)
                        DO 21 MU = 1,IMUE
                          CONE(MU,N) = CONE(MU,N)+X*S(MU)
   21                   CONTINUE
                      ENDIF
*
                    ENDIF
*
*   K = NU-1
*
                    IF (IBRE .EQ. 2) GOTO 25
*
                    IMUE = 4
                    IF (IBRE .EQ. 4) IMUE = 3
                    K = NU-1
*
                    IF (NCODE .EQ. N) THEN
                      X=XCODE
                    ELSE
                      ITKMO = NUD-2
                      IF (ITRIG (KS(1),KS(4),ITKMO) .EQ. 0) GOTO 23
                      IF (ITRIG (KS(2),KS(3),ITKMO) .EQ. 0) GOTO 23
                      J1(MJA) = ITKMO
                      CALL GENSUM (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,
     :                  JE9,KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                  JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     :                  KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,
     :                  INVE6J,X)
                    ENDIF
*
                    IF (IBUG2 .NE. 0) WRITE (99,306) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK (S,IS,KAPS,NU,K,IBRE,2)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,IMUE)
                      DO 22 MU = 1,IMUE
                        CONE(MU,N) = CONE(MU,N)+X*S(MU)
   22                 CONTINUE
                    ENDIF
*
                  ENDIF
*
*   K = NU+1
*
   23             IF ((IBRE .EQ. 2) .OR. (N .EQ. NE2)) GOTO 25
*
                  NCODE = N+1
                  XCODE = 0.0D 00
                  IMUE = 12
                  IF (IBRE .EQ. 4) IMUE = 7
                  ITKMO = NUD+2
*
                  IF ((ITRIG (KS(1),KS(4),ITKMO) .NE. 0) .AND.
     :                (ITRIG (KS(2),KS(3),ITKMO) .NE. 0)) THEN
                    K = NU+1
                    J1(MJA) = ITKMO
                      CALL GENSUM (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,
     :                  JE9,KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                  JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     :                  KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,
     :                  INVE6J,X)
                    XCODE = X
                    IF (IBUG2 .NE. 0) WRITE (99,306) NU,K,X
*
                    IF (ABS (X) .GE. EPS) THEN
                      X = X*PROD
                      CALL CXK (S,IS,KAPS,NU,K,IBRE,2)
                      IF (IBUG2 .NE. 0) WRITE (99,305)
     :                                  (S(III),III = 1,IMUE)
                      DO 24 MU = 1,IMUE
                        CONE (MU,N) = CONE(MU,N)+X*S(MU)
   24                 CONTINUE
                    ENDIF
*
                  ENDIF
*
   25           CONTINUE
*
              ENDIF
*
   26       CONTINUE
   27     CONTINUE
   28   CONTINUE
   29 CONTINUE
*
*   4.4 Insert outside factors
*
   30 IF (IBRD .GE. 0) THEN
        PRODD = KS(1)*KS(4)
        PRODD = CONST/SQRT (PRODD)
        IF ((IA1 .EQ. IB1) .AND. (IA2 .EQ. IB2))
     :     PRODD = 0.5D 00*PRODD
        DO 32 N = 1,ND2
          DO 31 MU = 1,12
            COND(MU,N) = COND(MU,N)*PRODD
   31     CONTINUE
   32   CONTINUE
      ENDIF
*
      IF (IBRE .GE. 0) THEN
        PRODE = KS(1)*KS(3)
        PRODE = -CONST/SQRT(PRODE)
        DO 34 N = 1,NE2
          DO 33 MU = 1,12
            CONE(MU,N) = CONE(MU,N)*PRODE
   33     CONTINUE
   34   CONTINUE
      ENDIF
*
*   5.0 Output results
*
      IF (IBRD .GE. 0) THEN
        DO 35 N = 1,ND2
          NU = ND1+2*(N-1)
          ITYPE = 1
          IF (IBRD .EQ. 2) ITYPE = 3
          IF (IBRD .EQ. 3) ITYPE = 4
          CALL TALK (JA,JB,NU,IA1,IA2,IB1,IB2,ITYPE,COND(1,N))
          IF (IBRD .GT. 1) GOTO 35
          CALL TALK (JA,JB,NU,IA2,IA1,IB2,IB1,ITYPE,COND(2,N))
          CALL TALK (JA,JB,NU,IA1,IA2,IB2,IB1,ITYPE,COND(3,N))
          CALL TALK (JA,JB,NU,IA2,IA1,IB1,IB2,ITYPE,COND(4,N))
          IF (N .EQ. ND2) GOTO 35
          NUP1 = NU+1
          ITYPE = 2
          CALL TALK (JA,JB,NUP1,IA1,IA2,IB1,IB2,ITYPE,COND(5,N))
          CALL TALK (JA,JB,NUP1,IB1,IB2,IA1,IA2,ITYPE,COND(6,N))
          CALL TALK (JA,JB,NUP1,IA2,IA1,IB2,IB1,ITYPE,COND(7,N))
          CALL TALK (JA,JB,NUP1,IB2,IB1,IA2,IA1,ITYPE,COND(8,N))
          CALL TALK (JA,JB,NUP1,IA1,IA2,IB2,IB1,ITYPE,COND(9,N))
          CALL TALK (JA,JB,NUP1,IB2,IB1,IA1,IA2,ITYPE,COND(10,N))
          CALL TALK (JA,JB,NUP1,IA2,IA1,IB1,IB2,ITYPE,COND(11,N))
          CALL TALK (JA,JB,NUP1,IB1,IB2,IA2,IA1,ITYPE,COND(12,N))
   35   CONTINUE
      ENDIF
*
      IF (IBRE .LT. 0) RETURN
*
      DO 36 N = 1,NE2
        NU = NE1+2*(N-1)
        IF (IBRE .NE. 4) THEN
          ITYPE = 1
          IF (IBRE .EQ. 2) ITYPE = 3
          CALL TALK (JA,JB,NU,IA1,IB2,IB1,IA2,ITYPE,CONE(1,N))
          IF (IBRE .EQ. 2) GOTO 36
          CALL TALK (JA,JB,NU,IB2,IA1,IA2,IB1,ITYPE,CONE(2,N))
          CALL TALK (JA,JB,NU,IA1,IB2,IA2,IB1,ITYPE,CONE(3,N))
          CALL TALK (JA,JB,NU,IB2,IA1,IB1,IA2,ITYPE,CONE(4,N))
          IF (N .EQ. NE2) GOTO 36
          NUP1 = NU+1
          ITYPE = 2
          CALL TALK (JA,JB,NUP1,IA1,IB2,IB1,IA2,ITYPE,CONE(5,N))
          CALL TALK (JA,JB,NUP1,IB1,IA2,IA1,IB2,ITYPE,CONE(6,N))
          CALL TALK (JA,JB,NUP1,IB2,IA1,IA2,IB1,ITYPE,CONE(7,N))
          CALL TALK (JA,JB,NUP1,IA2,IB1,IB2,IA1,ITYPE,CONE(8,N))
          CALL TALK (JA,JB,NUP1,IA1,IB2,IA2,IB1,ITYPE,CONE(9,N))
          CALL TALK (JA,JB,NUP1,IA2,IB1,IA1,IB2,ITYPE,CONE(10,N))
          CALL TALK (JA,JB,NUP1,IB2,IA1,IB1,IA2,ITYPE,CONE(11,N))
          CALL TALK (JA,JB,NUP1,IB1,IA2,IB2,IA1,ITYPE,CONE(12,N))
        ELSE
          ITYPE = 5
          CALL TALK (JA,JB,NU,IB1,IA1,IB1,IA1,ITYPE,CONE(1,N))
          CALL TALK (JA,JB,NU,IA1,IB1,IB1,IA1,ITYPE,CONE(2,N))
          CALL TALK (JA,JB,NU,IA1,IB1,IA1,IB1,ITYPE,CONE(3,N))
          IF (N .EQ. NE2) GOTO 36
          NUP1 = NU+1
          ITYPE = 6
          CALL TALK (JA,JB,NUP1,IA1,IB1,IA1,IB1,ITYPE,CONE(4,N))
          CALL TALK (JA,JB,NUP1,IB1,IA1,IB1,IA1,ITYPE,CONE(5,N))
          CALL TALK (JA,JB,NUP1,IA1,IB1,IB1,IA1,ITYPE,CONE(6,N))
          CALL TALK (JA,JB,NUP1,IB1,IA1,IA1,IB1,ITYPE,CONE(7,N))
        ENDIF
   36 CONTINUE
*
      RETURN
*
*   6.0 Fault diagnostic prints
*
   98 IF (IBUG2 .NE. 0) WRITE (99,300)
      RETURN
*
  300 FORMAT ('BREIT: Spectator quantum numbers not diagonal for',
     :        ' non-interacting shells')
  301 FORMAT ('BREIT: Increase second dimension of arrays',
     :        ' COND(MU,N) and CONE(MU,N) to the new value of NUMAX,'
     :       /' (at least ',1I3,').')
  302 FORMAT ('BREIT: Subshells ',4I5)
  303 FORMAT ('  ND1 ND2 NE1 NE2 IBRD IBRE',6I5)
  304 FORMAT ('  Direct NU K recoupling coef ',2I5,1P,D20.9)
  305 FORMAT (' S',1P,8D15.7)
  306 FORMAT ('  Exchange NU K recoupling coef ',2I5,1P,D20.9)
  307 FORMAT ('  Statistical factor ',1P,D20.9)
*
      END
