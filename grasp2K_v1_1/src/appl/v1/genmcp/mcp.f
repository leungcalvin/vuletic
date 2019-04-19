************************************************************************
*                                                                      *
      SUBROUTINE MCP (RESTRT)
*                                                                      *
*   This routine controls the computation  and storage of the values   *
*   and all indices of the angular coefficients                        *
*                                                                      *
*                                       k                              *
*                   T  (ab)            V  (abcd)                       *
*                    rs                 rs                             *
*                                                                      *
*   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
*   c and d are orbital sequence numbers.  r and s are configuration   *
*   state function indices.                                            *
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO,          *
*                        TNSRJJ.                                       *
*               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNJCUP,PNTJQS,PNTRIQ
      POINTER (PNJCUP,JCUPDUMMY)
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG,F0INT,LDBPA,LFORDR,LINCR,RESTRT
      CHARACTER*20 CNUM
      CHARACTER*2 CK,NH
*
      DIMENSION TSHELL(NNNW)
*
      POINTER (PLISTV,LLISTV(0:*))
*
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
*
      EXTERNAL COR,CORD
*
      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /DEBUGA/LDBPA(5)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DEFAULT/NDEF
     :      /FOPARM/ICCUT
     :      /MCPA/KMAX
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STAT/PNTJQS,PNJCUP
*
*   Set the encoding key and its square
*
CGG      PARAMETER (KEYORB =121)
      PARAMETER (KEY = KEYORB, KEYSQ = KEYORB*KEYORB)
*
*   Establish the cutoff criterion
*
      PARAMETER (CUTOFF = 1.0D-10)
*
      OPEN (20,FILE='sms.20',FORM='UNFORMATTED',
     :      STATUS='UNKNOWN')
*
*   Allocate storage to the array that stores the lengths
*   of the lists of V coefficients of each multipolarity
*
      CALL ALLOC (PLISTV,KMAX+1,4)
*
*   If this is a restart, determine the pair of CSFs with which
*   the computation should begin
*
      IF (RESTRT) THEN
*
         CALL FNDBEG (JASTRT,JBSTRT,INDEX,LLISTT,LLISTV)
*
      ELSE
*
         JASTRT = 1
         JBSTRT = 1
         INDEX = 0
*
*   Initialize the counters for the total number of T coefficients
*   and V coefficients of each multipolarity
*
         LLISTT = 0
         DO 1 K = 0,KMAX
            LLISTV(K) = 0
    1    CONTINUE
*
      ENDIF
*
*   Set the rank (zero) and parity (even) for the one-particle
*   coefficients
*
      KA = 0
      IOPAR = 1
*
*   INCOR is 0: only coefficients between open subshells required
*
      INCOR = 0
*
*   Allocate storage for the arrays in BUFFER
*
      CALL ALCBUF (1)
*
*   Write header to  .dbg  file if appropriate
*
      IF (LDBPA(2) .OR. LDBPA(3)) WRITE (99,300)
*
*   JA and JB respectively refer to the initial and final states
*   in the list of NCF configurations
*
      DO 5 JA = JASTRT,NCF
*
         IF ( DIAG .OR.
     :       (LFORDR .AND. (JA .GT. ICCUT)) ) THEN
            JBEND = JA
         ELSE
            JBEND = NCF
         ENDIF

*         write(*,*) 'Per',JA,JBSTRT,JBEND
*
         DO 4 JB = JBSTRT,JBEND
*
*   LINCR is .TRUE. if INDEX is to be incremented by 1; there
*   is always a diagonal element in each column
*
            IF (JB .NE. JA) THEN
               LINCR = .TRUE.
            ELSE
               INDEX = INDEX+1
               LINCR = .FALSE.
            ENDIF
*
            IF (JB .NE. JA) THEN
*
*   Call the MCT package to compute T coefficients
*
               CALL TNSRJJ (KA,IOPAR,JA,JB,IA,IB,TSHELL)
*
*   Write T coefficients that have magnitudes greater than the
*   cutoff criterion
*
               IF (IA .NE. 0) THEN
                  IF (IA .NE. IB) THEN
                     TCOEFF = TSHELL(1)
                     IF (ABS (TCOEFF) .GT. CUTOFF) THEN
                        IF (LDBPA(2)) WRITE (99,301) JB,JA,
     :                        NP(IA),NH(IA),NP(IB),NH(IB),TCOEFF
                        IF (LINCR) THEN
                           INDEX = INDEX+1
                           LINCR = .FALSE.
                        ENDIF
                        LLISTT = LLISTT+1
                        LAB = MIN (IA,IB)*KEY+MAX (IA,IB)
                        WRITE (31) JA,INDEX,LAB,TCOEFF
                     ENDIF
                  ENDIF
               ENDIF
*
            ENDIF
*
*   Initialize
*
            NVCOEF = 0
*
*   Call the MCP package to generate V coefficients; ac and bd
*   are the density pairs
*
            CALL RKCO (JA,JB,COR,CORD,INCOR)
*
            DO 2 I = 1,NVCOEF
               VCOEFF = COEFF(I)
               IF (ABS (VCOEFF) .GT. CUTOFF) THEN
                  IA = LABEL(1,I)
                  IB = LABEL(2,I)
                  IC = LABEL(3,I)
                  ID = LABEL(4,I)
                  K  = LABEL(5,I)
                  F0INT = (K .EQ. 0) .AND.
     :                    (IA .EQ. IC) .AND.
     :                    (IB .EQ. ID)
                  IF (.NOT. F0INT) THEN
                     IF (LDBPA(3)) WRITE (99,302) K,JB,JA,
     :                  NP(IA),NH(IA),NP(IB),NH(IB),
     :                  NP(IC),NH(IC),NP(ID),NH(ID),VCOEFF
                     NSWAP = 0
                     IF (IA .GT. IC) THEN
                        ISWAP = IC
                        IC = IA
                        IA = ISWAP
                        NSWAP = NSWAP + 1
                     ENDIF
                     IF (IB .GT. ID) THEN
                        ISWAP = ID
                        ID = IB
                        IB = ISWAP
                        NSWAP = NSWAP + 1
                     ENDIF
                     IF (LINCR) THEN
                        INDEX = INDEX+1
                        LINCR = .FALSE.
                     ENDIF
                     LLISTV(K) = LLISTV(K)+1
                     LAC = IA*KEY+IC
                     LBD = IB*KEY+ID
                     IF (LAC .LT. LBD) THEN
                        LAB = LAC*KEYSQ+LBD
                     ELSE
                        LAB = LBD*KEYSQ+LAC
                     ENDIF
                     WRITE (32+K) JA,INDEX,LAB,VCOEFF
                     IF (K.EQ.1) WRITE (20) NSWAP
                  ENDIF
               ENDIF
    2       CONTINUE
*
*   All angular coefficients for this pair of CSFs have been
*   generated; update file 30
*
            IF ((JB .EQ. JA) .OR. (.NOT. LINCR))
     :         WRITE (30) JA,JB,INDEX
*
    4    CONTINUE
*
         CALL CONVRT (JA,CNUM,LCNUM)
         PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
*
         JBSTRT = JA+1
*
    5 CONTINUE
*
*   Deallocate storage that is no longer required
*
      CALL DALLOC (PNTRIQ)
      CALL DALLOC (PNTJQS)
      CALL DALLOC (PNJCUP)
      CALL ALCBUF (3)
*
*   Write out a report for this run
*
      IF (NDEF.NE.0) THEN
         WRITE (24,*)
         CALL CONVRT (LLISTT,CNUM,LENTH)
         WRITE (24,*) CNUM(1:LENTH)//' T coefficients generated;'
         DO 6 K = 0,KMAX
            CALL CONVRT (LLISTV(K),CNUM,LENTH)
            CALL CONVRT (K,CK,LCK)
            WRITE (24,*) CNUM(1:LENTH)//' V(k='//CK(1:LCK)
     :      //') coefficients generated;'
    6    CONTINUE
      ENDIF
*
*   Set up sparse structure definition arrays in file 30
*
      CALL SETSDA (INDEX,LDBPA(4))
*
*   Sort MCP coefficients into integral-based lists
*
      IF (NDEF.NE.0) WRITE (24,*)
      DO 7 I = 31,32+KMAX
         IF (I .GT. 31) THEN
            CALL SORT (I,LLISTV(I-32),NTGI,LDBPA(3))
            CALL CONVRT (I-32,CK,LCK)
            CALL CONVRT (NTGI,CNUM,LENTH)
            IF (NDEF.NE.0)
     :         WRITE (24,*) ' k = '//CK(1:LCK)//': '
     :         //CNUM(1:LENTH)//' Slater integrals;'
         ELSE
            CALL SORT (I,LLISTT,NTGI,LDBPA(2))
            CALL CONVRT (NTGI,CNUM,LENTH)
            IF (NDEF.NE.0)
     :      WRITE (24,*) CNUM(1:LENTH)//' I(ab) integrals;'
         ENDIF
    7 CONTINUE
      CALL DALLOC (PLISTV)
*
      RETURN
*
  300 FORMAT (/'From MCP:')
  301 FORMAT (' T_[',1I3,',',1I3,']',
     :   ' (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
  302 FORMAT (' V^[(',1I2,')]_[',1I3,',',1I3,']',
     :   ' (',1I2,1A2,',',1I2,1A2,';',
     :        1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
*
      END
