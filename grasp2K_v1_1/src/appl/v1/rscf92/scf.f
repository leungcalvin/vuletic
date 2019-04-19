************************************************************************
*                                                                      *
      SUBROUTINE SCF (EOL)
*                                                                      *
*   This  subroutine  performs  the SCF iterations. The procedure is   *
*   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
*   (1986) 290.                                                        *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
*               [RSCF92]: IMPROV, MATRIX, MAXARR, NEWCO, ORBOUT,       *
*                         ORTHSC, SETLAG.                              *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last update: 22 Dec 1992    *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PCCMIN,CCMINDUMMY),
     :        (PNEVAL,EVALDUMMY),
     :        (PNEVEC,EVECDUMMY),
     :        (PNTEMT,EMTDUMMY),
     :        (PNIECC,IECCDUMMY), (PNTECV,ECVDUMMY),
     :        (PNTRIQ,RIQDUMMY), 
     :        (PNTNDA,NDADUMMY), (PNTNXA,NXADUMMY),
     :        (PNTNYA,NYADUMMY), (PNTRDA,RDADUMMY),
     :        (PNTRXA,RXADUMMY), (PNTRYA,RYADUMMY),
     :        (PIATJP,IATJPDUMMY), (PIASPA,IASPADUMMY),
     :        (PNTJQS,JQSDUMMY), (PNJCUP,JCUPDUMMY)
      LOGICAL CONVG,EOL,LDBPR,LFIX,orthst,lsort
*
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /LAGR/PNTECV,PNIECC,NEC
     :      /MCPA/KMAXF
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORBA/IORDER(NNNW)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /SYMA/PIATJP,PIASPA
     :      /STAT/PNTJQS,PNJCUP
     :      /ORTHCT/orthst

      COMMON/DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde

      IF (ndef .EQ. 0) THEN
         lsort = .FALSE.
      ELSE
  123    WRITE(istde,*) 'Orthonomalization order? '
         WRITE(istde,*) '     1 -- Update order'
         WRITE(istde,*) '     2 -- Self consistency connected'
         READ(istdi,*) j
         IF (j .EQ. 1) THEN
            lsort = .FALSE.
         ELSE IF( j .EQ. 2 ) THEN
            lsort = .TRUE.
         ELSE
            WRITE(istde,*) 'Input is wrong, redo...'
            GOTO 123
         ENDIF
      END IF
*
*   Deallocate storage that will no longer be used
*
      CALL DALLOC (PNTJQS)
*
*   Initialisations for (E)OL calculations
*
      IF (EOL) THEN
         CALL ALLOC (PNEVAL,NCMIN,8)
         CALL ALLOC (PNEVEC,NCF*NCMIN,8)
         CALL ALLOC (PIATJP,NCMIN,4)
         CALL ALLOC (PIASPA,NCMIN,4)
      ENDIF
*
*   General initializations
*
      READ (30) NELMNT
      CALL ALLOC (PIENDC,NCF+1,4)
      CALL ALLOC (PNIROW,NELMNT,4)
      READ (30) (IENDC(I),I = 0,NCF),(IROW(I),I = 1,NELMNT)
      CLOSE (30)
*
      NDDIM = 0
      NXDIM = 0
      NYDIM = 0
*
      CALL SETLAG (EOL)
*
*
*   For (E)OL calculations, determine the level energies and
*   mixing coefficients
*
      IF (EOL) THEN
         CALL MATRIX
         CALL NEWCO
      ENDIF

      DO 4 NIT = 1,NSCF
!always         WRITE(istde,*) ' Iteration number ', NIT
         WRITE (*,301) NIT
*
*   For all pairs constrained through a Lagrange multiplier, compute
*   the Lagrange multiplier
*
         CALL SETLAG (EOL)
*
*   Improve all orbitals in turn
*
         WRITE (*,302) 
         DO 1 J = 1,NW
            JSEQ = IORDER(J)
            IF (.NOT. LFIX(JSEQ)) THEN
               CALL IMPROV (EOL,JSEQ,lsort)
            ENDIF
    1    CONTINUE
*
*   For KOUNT = 1 to NSIC: find the least self-consistent orbital;
*   improve it
*
!         write(istde,*) 'nsic=',nsic
         DO 2 KOUNT = 1,NSIC
            CALL MAXARR (K)
            IF (K .EQ. 0) THEN
               CONVG = .TRUE.
               GOTO 3
            ELSE
               IF (SCNSTY(K) .LE. ACCY) THEN
                  CONVG = .TRUE.
                  GOTO 3
               ENDIF
            ENDIF
            CALL IMPROV (EOL,K,lsort)
    2    CONTINUE
*
!XHH         CALL ORBOUT
*
         CALL MAXARR (K)
*
         IF (K .EQ. 0) THEN
            CONVG = .TRUE.
         ELSE
            IF (SCNSTY(K) .LE. ACCY) THEN
               CONVG = .TRUE.
            ELSE
               CONVG = .FALSE.
            ENDIF
         ENDIF
*
    3    IF (LDBPR(24)) CALL PRWF (0)
*
*   Perform Gram-Schmidt process
*
!   For OL calculation, orthst is true and orbitals are orthonormalized
!   in subroutine improv. For AL calculation, orthst is false.
         IF( .NOT. orthst ) CALL ORTHSC
*
*   Write the subshell radial wavefunctions to the .rwf file
*
         CALL ORBOUT
*
         IF (EOL) THEN
            CALL MATRIX
            CALL NEWCO
         ENDIF
         IF (CONVG) THEN
            IF (LDBPR(25) .AND. (.NOT. LDBPR(24))) CALL PRWF (0)
!            IF (EOL) CALL MATRIX
            GOTO 5
         ENDIF
*
    4 CONTINUE
*
      WRITE(istde,*) ' Maximum iterations in SCF Exceeded.'
*
*   Close MCP coefficient files
*
    5 DO 6 I = 31,32+KMAXF
         CLOSE (I)
    6 CONTINUE
*
*   Close the .rwf file
*
      CLOSE (23)
*
*   Close the .mix file
*
      CLOSE (25)
*
*   Deallocate storage
*
      IF (NEC .GT. 0) THEN
         CALL DALLOC (PNIECC)
         CALL DALLOC (PNTECV)
         CALL DALLOC (PNTRIQ)
      ENDIF
      IF (NDDIM .GT. 0) THEN
         CALL DALLOC (PNTRDA)
         CALL DALLOC (PNTNDA)
      ENDIF
      IF (NXDIM .GT. 0) THEN
         CALL DALLOC (PNTRXA)
         CALL DALLOC (PNTNXA)
      ENDIF
      IF (NYDIM .GT. 0) THEN
         CALL DALLOC (PNTRYA)
         CALL DALLOC (PNTNYA)
      ENDIF
      IF (EOL) THEN
         CALL DALLOC (PIENDC)
         CALL DALLOC (PNIROW)
      ENDIF
*
      RETURN
*
  301 FORMAT (/' Iteration number ',1I3
     :        /' ----------------')
  302 FORMAT (74X,'Self-     Damping'
     :        /'Subshell      Energy              P0         '
     :        ,'      Norm        Method  consistency   factor    JP'
     :        ,'   MTP   INV   NNP'/)
*
      END
