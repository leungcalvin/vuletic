************************************************************************
*                                                                      *
      SUBROUTINE SETLAG (EOL)
*                                                                      *
*   Sets up the data structure  pertaining to the Lagrange multipli-   *
*   ers  on the first entry;  on subsequent calls it  determines new   *
*   estimates for the multipliers.                                     *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, QUAD, RINTI,                           *
*               [RSCF92]: DACON, SETCOF, XPOT,  YPOT.                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
      COMMON/CORE/NCORE
!$Id: setlag.f,v 1.1 2003/09/30 05:51:40 georgio Exp $
!$Log: setlag.f,v $
!Revision 1.1  2003/09/30 05:51:40  georgio
!
!added
!
!Revision 1.3  1997/06/02 21:58:33  xhh
!*** empty log message ***
!
!Revision 1.2  1997/03/04 00:37:08  xhh
!Lagrange multipliers (see NCORE)
!
! Changes made to keep Lagrange multipliers among closed-shell electrons
! zero (suggested by CFF):
! 1) Added a common statement COMMON/CORE/NCORE to read the value NCORE
!   obtained through RSCF92(COMMON) <-- SETCSL(Dummy) <-- LODCSL(Dummy)
! 2) Replaced  LIP1 = LI+1  by  LIP1 = MAX(NCORE,LI)+1
! XHH 1997.01.23
! Applied the order of updating orbitals to the Lagrange multilier
! calculation (look iorder() for changes)
! 97.03.14

Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*2 NH
      LOGICAL EOL,FIRST,FIXLI,FIXLJ,FULLI,FULLJ,LFIX
*
      POINTER (PNTECV,ECV(1))
      POINTER (PNIECC,IECC(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      DIMENSION YPJ(NNNP),YPM(NNNP),
     :          XPJ(NNNP),XPM(NNNP),
     :          XQJ(NNNP),XQM(NNNP)
*
      COMMON/DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /LAGR/PNTECV,PNIECC,NEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF1/UCF(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /ORBA/IORDER(NNNW)
*
      DATA FIRST/.TRUE./
*
Cww      PARAMETER (P001 = 1.0D-03)
      PARAMETER (P001 = 1.0D-01)
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
Cww      EPS = ACCY*0.1D 00
      EPS = ACCY*0.01D 00
*
      IF (FIRST) THEN
*
*   Initialize
*
         NEC = 0
*
*   Determine the total number of Lagrange multipliers
*
         NLMDIM = 0
*
         NWM1 = NW-1
!         DO 2 LI = 1,NWM1
         DO 2 LIraw = 1, NWM1
            LI = iorder(LIraw)
            !LIP1 = LI+1   ! Changed to the next one
            LIP1 = MAX(NCORE,LIraw)+1
            NAKLI = NAK(LI)
            FIXLI = LFIX(LI)
            FULLI = ABS ( UCF(LI)-DBLE (NKJ(LI)+1) )
     :              .LT. EPS
!            DO 1 LJ = LIP1,NW
            DO 1 LJraw = LIP1, NW
               LJ = iorder(LJraw)
               FIXLJ = LFIX(LJ)
               FULLJ = ABS ( UCF(LJ)-DBLE (NKJ(LJ)+1) )
     :                 .LT. EPS
               IF ( (NAK(LJ) .EQ. NAKLI) .AND.
     :              (.NOT. (FIXLI .AND. FIXLJ)) .AND.
     :              (.NOT. (FULLI .AND. FULLJ)) ) THEN
                  NEC = NEC+1
                  IF (NEC .GT. NLMDIM) THEN
                     IF (NLMDIM .EQ. 0) THEN
                        CALL ALLOC (PNTECV,NEC,8)
                        CALL ALLOC (PNIECC,NEC,4)
                     ELSE
                        CALL RALLOC (PNTECV,NLMDIM,NEC,8)
                        CALL RALLOC (PNIECC,NLMDIM,NEC,4)
                     ENDIF
                     NLMDIM = NEC
                  ENDIF
*
*   Encode index
*
                  IECC(NEC) = LI+KEY*LJ
               ENDIF
    1       CONTINUE
    2    CONTINUE
*
*   Print information about Lagrange multipliers
*
         IF (NEC .EQ. 0) THEN
            WRITE (*,302)
         ELSE
            WRITE (*,304)
            DO 3 LI = 1,NEC
*
*   Decode index
*
               IECCLI = IECC(LI)
               L1 = IECCLI/KEY
               L2 = IECCLI-KEY*L1
*
               WRITE (*,305) NP(L2),NH(L2),NP(L1),NH(L1)
    3       CONTINUE
         ENDIF
*
         FIRST = .FALSE.
*
         RETURN
*
      ENDIF
*
*   Return if there are no Lagrange multipliers to be computed
*
      IF (NEC .EQ. 0) RETURN
*
      WRITE (*,306)
*
*   Compute Lagrange multipliers
*
      JLAST = 0
      MLAST = 0
*
      DO 10 LI = 1,NEC
*
*   Decode index
*
         IECCLI = IECC(LI)
         M = IECCLI/KEY
         J = IECCLI-KEY*M
*
         IF (J .NE. JLAST) THEN
            UCFJ = UCF(J)
            CALL SETCOF (EOL,J)
            CALL YPOT (J)
            CALL XPOT (J)
            CALL DACON
            DO 4 I = 1,N
               YPJ(I) = YP(I)
               XPJ(I) = XP(I)
               XQJ(I) = XQ(I)
    4       CONTINUE
            JLAST = J
         ENDIF
*
         IF (M .NE. MLAST) THEN
            UCFM = UCF(M)
            CALL SETCOF (EOL,M)
            CALL YPOT (M)
            CALL XPOT (M)
            CALL DACON
            DO 5 I = 1,N
               YPM(I) = YP(I)
               XPM(I) = XP(I)
               XQM(I) = XQ(I)
    5       CONTINUE
            MLAST = M
         ENDIF
*
         MTP = MAX (MF(J),MF(M))
*
         IF (LFIX(M)) THEN
            TA(1) = 0.0D 00
            DO 6 I = 2,MTP
               TA(I) = RPOR(I)*( ( PF(I,M)*XQJ(I)
     :                            -QF(I,M)*XPJ(I) )
     :                          *C
     :                          +( PF(I,M)*PF(I,J)
     :                            +QF(I,M)*QF(I,J) ) * YPJ(I) )
    6       CONTINUE
            CALL QUAD (RESULT)
            RIMJ = RINTI (M,J,1)
            ECV(LI) = (RESULT-RIMJ)*UCFJ
         ELSEIF (LFIX(J)) THEN
            TA(1) = 0.0D 00
            DO 7 I = 2,MTP
               TA(I) = RPOR(I)*( ( PF(I,J)*XQM(I)
     :                            -QF(I,J)*XPM(I) )
     :                          *C
     :                          +( PF(I,J)*PF(I,M)
     :                            +QF(I,J)*QF(I,M) ) * YPM(I) )
    7       CONTINUE
            CALL QUAD (RESULT)
            RIJM = RINTI (J,M,1)
            ECV(LI) = (RESULT-RIJM)*UCFM
         ELSE
            QDIF = ABS ((UCFJ-UCFM)/MAX (UCFJ,UCFM))
            IF (QDIF .GT. P001) THEN
               OBQDIF = 1.0D 00/UCFJ-1.0D 00/UCFM
               TA(1) = 0.0D 00
               DO 8 I = 2,MTP
                  TA(I) = RPOR(I)*( ( PF(I,M)*XQJ(I)
     :                               -QF(I,M)*XPJ(I)
     :                               -PF(I,J)*XQM(I)
     :                               +QF(I,J)*XPM(I) )
     :                             *C
     :                             +( YPJ(I)-YPM(I) )
     :                             *( PF(I,M)*PF(I,J)
     :                               +QF(I,M)*QF(I,J) ) )
    8          CONTINUE
               CALL QUAD (RESULT)
               ECV(LI) = RESULT/OBQDIF
            ELSE
               OBQSUM = 1.0D 00/UCFJ+1.0D 00/UCFM
               TA(1) = 0.0D 00
               DO 9 I = 2,MTP
                  TA(I) = RPOR(I)*( ( PF(I,M)*XQJ(I)
     :                               -QF(I,M)*XPJ(I)
     :                               +PF(I,J)*XQM(I)
     :                               -QF(I,J)*XPM(I) )
     :                             *C
     :                             +( YPJ(I)+YPM(I) )
     :                             *( PF(I,M)*PF(I,J)
     :                               +QF(I,M)*QF(I,J) ) )
    9          CONTINUE
               CALL QUAD (RESULT)
               RIMJ = RINTI (M,J,1)
               ECV(LI) = (RESULT-2.0D 00*RIMJ)/OBQSUM
            ENDIF
         ENDIF
*
         WRITE (*,307) NP(J),NH(J),NP(M),NH(M),ECV(LI)
*
   10 CONTINUE
*
      RETURN
*
  302 FORMAT (/'Lagrange multipliers are not required')
  304 FORMAT (/'Include Lagrange multipliers between:'/)
  305 FORMAT (13X,2(2X,1I2,1A2))
  306 FORMAT (/'Lagrange multipliers:'/)
  307 FORMAT (13X,2(2X,1I2,1A2),2X,1PD16.9)
*
      END
