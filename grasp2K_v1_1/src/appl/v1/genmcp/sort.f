************************************************************************
*                                                                      *
      SUBROUTINE SORT (NFILE,NCOEFF,NTGRAL,LPRINT)
*                                                                      *
*   This routine sorts lists                                           *
*                                                                      *
*                     (ICLMN,INDEX,LABEL,COEFF)                        *
*   into lists                                                         *
*                     (LABEL,ICLMN,INDEX,COEFF)                        *
*                                                                      *
*   using Heapsort. File NFILE is closed by this routine.  NCOEFF is   *
*   the number of triads (INDEX, ...) and is an input. NTGRAL is the   *
*   number of different values of LABEL, and is an output. If LPRINT   *
*   is  .TRUE. , the contents of the sorted file are interpreted and   *
*   printed to unit 99.                                                *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL DIAG,LFORDR,LPRINT
      CHARACTER*20 CNUM
      CHARACTER*8 SRTLAB
      CHARACTER*3 MCPLAB
      CHARACTER*2 CK,NH
*
      POINTER (PCOEFF,COEFF(1)), (PICLMN,ICLMN(1)), (PINDEX,INDEX(1)),
     : (PLABEL,LABEL(1)), (PNSWAP,NSWAP(1))
*
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
*
*   This is the encoding key
*
CGG      PARAMETER (KEYORB =121)
      PARAMETER (KEY = KEYORB)
*
      REWIND (NFILE)
      IF (NFILE.EQ.33) REWIND (20)
*
*   Has the file already been sorted? Is it empty other than
*   the header? Skip the sort if so.
*
      READ (NFILE) MCPLAB,SRTLAB
      IF (SRTLAB .EQ. '  SORTED') GOTO 6
*
*   Message
*
      CALL CONVRT (NCOEFF,CNUM,LCNUM)
      IF (NFILE .GT. 31) THEN
         CALL CONVRT (NFILE-32,CK,LCK)
         PRINT *, 'Sorting '//CNUM(1:LCNUM)
     :       //' V(k='//CK(1:LCK)//') coefficients ...'
      ELSE
         PRINT *, 'Sorting '//CNUM(1:LCNUM)//' T coefficients ...'
      ENDIF
*
*   Position the file
*
      READ (NFILE) NELEC,NCF,NW
      READ (NFILE) DIAG,ICCUT,LFORDR
*
*   Sort the list
*
      IF (NCOEFF .GT. 0) THEN
*
*   Allocate storage for all required arrays
*
         CALL ALLOC (PCOEFF,NCOEFF,8)
         CALL ALLOC (PICLMN,NCOEFF,4)
         CALL ALLOC (PINDEX,NCOEFF,4)
         CALL ALLOC (PLABEL,NCOEFF,4)
         IF (NFILE.EQ.33) CALL ALLOC (PNSWAP,NCOEFF,4)
*
*   Read arrays into memory from NFILE
*
         DO 1 I = 1,NCOEFF
            READ (NFILE) ICLMN(I),INDEX(I),LABEL(I),COEFF(I)
    1    CONTINUE

         IF (NFILE.EQ.33) THEN
            DO 11 I = 1,NCOEFF
               READ (20) NSWAP(I)
   11       CONTINUE
         ENDIF                
*
      ENDIF
*
*   Sort LABEL into ascending order using the heapsort algorithm;
*   move the associated members of COEFF and INDEX in the same
*   manner; the code below is adapted from Press et al.
*
      IF (NFILE.EQ.33) THEN
      IF (NCOEFF .GT. 1) THEN
*
         L = NCOEFF/2+1
         IR = NCOEFF
    2    IF (L .GT. 1) THEN
            L = L-1
            COF = COEFF(L)
            ICL = ICLMN(L)
            IND = INDEX(L)
            LAB = LABEL(L)
            NSW = NSWAP(L)
         ELSE
            COF = COEFF(IR)
            ICL = ICLMN(IR)
            IND = INDEX(IR)
            LAB = LABEL(IR)
            NSW = NSWAP(IR)
            COEFF(IR) = COEFF(1)
            ICLMN(IR) = ICLMN(1)
            INDEX(IR) = INDEX(1)
            LABEL(IR) = LABEL(1)
            NSWAP(IR) = NSWAP(1)
            IR = IR-1
            IF (IR .EQ. 1) THEN
               COEFF(1) = COF
               ICLMN(1) = ICL
               INDEX(1) = IND
               LABEL(1) = LAB
               NSWAP(1) = NSW
               GOTO 4
            ENDIF
         ENDIF
         I = L
         J = L+L
    3    IF (J .LE. IR) THEN
            IF (J .LT. IR) THEN
               IF (LABEL(J) .LT. LABEL(J+1)) J = J+1
            ENDIF
            IF (LAB .LT. LABEL(J)) THEN
               COEFF(I) = COEFF(J)
               ICLMN(I) = ICLMN(J)
               INDEX(I) = INDEX(J)
               LABEL(I) = LABEL(J)
               NSWAP(I) = NSWAP(J)
               I = J
               J = J+J
            ELSE
               J = IR+1
            ENDIF
            GOTO 3
         ENDIF
         COEFF(I) = COF
         ICLMN(I) = ICL
         INDEX(I) = IND
         LABEL(I) = LAB
         NSWAP(I) = NSW
         GOTO 2
*
      ENDIF
      ELSE
*
*   Sort LABEL into ascending order using the heapsort algorithm;
*   move the associated members of COEFF and INDEX in the same
*   manner; the code below is adapted from Press et al.
*
      IF (NCOEFF .GT. 1) THEN
*
         L = NCOEFF/2+1
         IR = NCOEFF
   92    IF (L .GT. 1) THEN
            L = L-1
            COF = COEFF(L)
            ICL = ICLMN(L)
            IND = INDEX(L)
            LAB = LABEL(L)
         ELSE
            COF = COEFF(IR)
            ICL = ICLMN(IR)
            IND = INDEX(IR)
            LAB = LABEL(IR)
            COEFF(IR) = COEFF(1)
            ICLMN(IR) = ICLMN(1)
            INDEX(IR) = INDEX(1)
            LABEL(IR) = LABEL(1)
            IR = IR-1
            IF (IR .EQ. 1) THEN
               COEFF(1) = COF
               ICLMN(1) = ICL
               INDEX(1) = IND
               LABEL(1) = LAB
               GOTO 4
            ENDIF
         ENDIF
         I = L
         J = L+L
   93    IF (J .LE. IR) THEN
            IF (J .LT. IR) THEN
               IF (LABEL(J) .LT. LABEL(J+1)) J = J+1
            ENDIF
            IF (LAB .LT. LABEL(J)) THEN
               COEFF(I) = COEFF(J)
               ICLMN(I) = ICLMN(J)
               INDEX(I) = INDEX(J)
               LABEL(I) = LABEL(J)
               I = J
               J = J+J
            ELSE
               J = IR+1
            ENDIF
            GOTO 93
         ENDIF
         COEFF(I) = COF
         ICLMN(I) = ICL
         INDEX(I) = IND
         LABEL(I) = LAB
         GOTO 92
*
      ENDIF
      ENDIF
*
*   Sorting complete; rewrite the file header
*
    4 REWIND (NFILE)
      IF (NFILE.EQ.33) REWIND (20)
      WRITE (NFILE) 'MCP','  SORTED'
      WRITE (NFILE) NELEC,NCF,NW
      WRITE (NFILE) DIAG,ICCUT,LFORDR
*
*   Write the sorted list to NFILE
*
      IF (NCOEFF .GT. 0) THEN
*
         LAST = LABEL(1)
         IBEG = 1
         IEND = 1
         NTGRAL = 1
*
         DO 5 I = 2,NCOEFF
            IF (LABEL(I) .EQ. LAST) THEN
               IEND = IEND+1
            ELSE
               WRITE (NFILE+40,*) LAST,IEND-IBEG+1
               WRITE (NFILE) LAST,IEND-IBEG+1
               WRITE (NFILE)
     :            (ICLMN(J),INDEX(J),COEFF(J),J = IBEG,IEND)
               IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,IEND)
               NTGRAL = NTGRAL+1
               LAST = LABEL(I)
               IBEG = IEND+1
               IEND = IBEG
            ENDIF
    5    CONTINUE
*
         IF (IBEG .LE. NCOEFF) THEN
            WRITE (NFILE+40,*) LAST,NCOEFF-IBEG+1
            WRITE (NFILE) LAST,NCOEFF-IBEG+1
            WRITE (NFILE) (ICLMN(J),INDEX(J),COEFF(J),J = IBEG,NCOEFF)
            IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,NCOEFF)
         ENDIF
*
*   Deallocate storage
*
         CALL DALLOC (PCOEFF)
         CALL DALLOC (PICLMN)
         CALL DALLOC (PINDEX)
         CALL DALLOC (PLABEL)
         IF (NFILE.EQ.33) CALL DALLOC (PNSWAP)
*
      ELSE
*
         NTGRAL = 0
*
      ENDIF
*
*   Completion message
*
      CALL CONVRT (NTGRAL,CNUM,LCNUM)
      PRINT *, ' ... sort complete; '//CNUM(1:LCNUM)//' integrals;'
*
*   Debug printout
*
    6 IF (LPRINT) THEN
         WRITE (99,300)
         NDIM = 1
         CALL ALLOC (PCOEFF,NDIM,8)
         CALL ALLOC (PICLMN,NDIM,4)
         CALL ALLOC (PINDEX,NDIM,4)
         REWIND (NFILE)
         READ (NFILE)
         READ (NFILE)
         READ (NFILE)
         NTGRAL = 0
         IF (NFILE .EQ. 31) THEN
    7       READ (NFILE,IOSTAT = IOS) LAB,NCONTR
            IF (IOS .EQ. 0) THEN
               NTGRAL = NTGRAL+1
               IA = MOD (LAB,KEY)
               IB = LAB/KEY
               WRITE (99,301) NP(IA),NH(IA),NP(IB),NH(IB)
               IF (NCONTR .GT. NDIM) THEN
                  CALL DALLOC (PCOEFF)
                  CALL DALLOC (PICLMN)
                  CALL DALLOC (PINDEX)
                  NDIM = NCONTR
                  CALL ALLOC (PCOEFF,NDIM,8)
                  CALL ALLOC (PICLMN,NDIM,4)
                  CALL ALLOC (PINDEX,NDIM,4)
               ENDIF
               READ (NFILE) (ICLMN(J),INDEX(J),COEFF(J),J = 1,NCONTR)
               DO 8 J = 1,NCONTR
                  WRITE (99,302) ICLMN(J),INDEX(J),COEFF(J)
    8          CONTINUE
               GOTO 7
            ENDIF
            WRITE (99,303) NTGRAL
         ELSE
            K = NFILE-32
    9       READ (NFILE,IOSTAT = IOS) LAB,NCONTR
            IF (IOS .EQ. 0) THEN
               NTGRAL = NTGRAL+1
               ID = MOD (LAB,KEY)
               LAB = LAB/KEY
               IB = MOD (LAB,KEY)
               LAB = LAB/KEY
               IC = MOD (LAB,KEY)
               IA = LAB/KEY
               WRITE (99,304) K,NP(IA),NH(IA),NP(IB),NH(IB),
     :                          NP(IC),NH(IC),NP(ID),NH(ID)
               IF (NCONTR .GT. NDIM) THEN
                  CALL DALLOC (PCOEFF)
                  CALL DALLOC (PICLMN)
                  CALL DALLOC (PINDEX)
                  NDIM = NCONTR
                  CALL ALLOC (PCOEFF,NDIM,8)
                  CALL ALLOC (PICLMN,NDIM,4)
                  CALL ALLOC (PINDEX,NDIM,4)
               ENDIF
               READ (NFILE) (ICLMN(J),INDEX(J),COEFF(J),J = 1,NCONTR)
               DO 10 J = 1,NCONTR
                  WRITE (99,305) K,ICLMN(J),INDEX(J),COEFF(J)
   10          CONTINUE
               GOTO 9
            ENDIF
            WRITE (99,303) NTGRAL
         ENDIF
      ENDIF
*
*   Close NFILE; it is no longer required by GENMCP
*
      CLOSE (NFILE)
      IF (NFILE.EQ.33) CLOSE (20)
*
      RETURN
*
  300 FORMAT (/'From SORT:')
  301 FORMAT (' I(',1I2,1A2,',',1I2,1A2,'):')
  302 FORMAT ('  T_[',1I2,',',1I4,'] = ',1PD19.12)
  303 FORMAT ('  Number of integrals is ',1I4)
  304 FORMAT (' R^[(',1I2,')] (',1I2,1A2,',',1I2,1A2,';'
     :                          ,1I2,1A2,',',1I2,1A2,'):')
  305 FORMAT ('  V^[(',1I2,')]_[',1I2,',',1I4,'] = ',1PD19.12)
*
      END
