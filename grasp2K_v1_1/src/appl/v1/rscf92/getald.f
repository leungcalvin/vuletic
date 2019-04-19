************************************************************************
*                                                                      *
      SUBROUTINE GETALD
*                                                                      *
*   Interactively determines the data governing AL problem.            *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, IQ, ITJPO, GETYN.                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PCCMIN,PWEIGH,PNTRIQ
      POINTER (PCCMIN,CCMINDUMMY)
      POINTER (PWEIGH,WEIGHDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL GETYN,LFIX,ORTHST,YES
*
      POINTER (PNTRWT,WT(1))
*
      COMMON/DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORTHCT/ORTHST
     :      /SCF1/UCF(NNNW)
*
      NCMIN = 0
*
      CALL ALLOC (PNTRWT,NCF,8)
*
      PRINT *, 'Standard ASF weights are (2J+1)'
      PRINT *, ' values; revise?'
      YES = GETYN ()
*
      IF (YES) THEN
         PRINT *, 'Assign all ASFs the same weight?'
         YES = GETYN ()
         IF (YES) THEN
            DO 1 I = 1,NCF
               WT(I) = 1.0D 00
    1       CONTINUE
            SUM = DBLE (NCF)
         ELSE
    2       PRINT *, 'Enter the (relative) weights for'
            PRINT *, ' all ASFs:'
            READ (*,*) (WT(I),I = 1,NCF)
            SUM = 0.0D 00
            DO 3 I = 1,NCF
               IF (WT(I) .LE. 0.0D 00) THEN
                  PRINT *, 'GETALD: Weights must exceed 0;'
                  GOTO 2
               ELSE
                  SUM = SUM+WT(I)
               ENDIF
    3       CONTINUE
         ENDIF
      ELSE
         SUM = 0.0D 00
         DO 4 I = 1,NCF
            FTJPOI = DBLE (ITJPO (I))
            WT(I) = FTJPOI
            SUM = SUM+FTJPOI
    4    CONTINUE
      ENDIF
*
      SUM = 1.0D 00/SUM
      DO 5 I = 1,NCF
         WT(I) = SUM*WT(I)
    5 CONTINUE
*
      DO 7 J = 1,NW
         SUM = 0.0D 00
         DO 6 I = 1,NCF
            SUM = SUM+WT(I)*DBLE (IQ (J,I))
    6    CONTINUE
         UCF(J) = SUM
    7 CONTINUE
*
      NSCF = 12
      NSIC = 2+(NW-NFIX)/4
      ORTHST = .FALSE.
*
      RETURN
      END
