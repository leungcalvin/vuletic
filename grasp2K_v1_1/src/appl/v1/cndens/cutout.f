************************************************************************
*                                                                      *
      SUBROUTINE CUTOUT(NAME)
*                                                                      *
*   Condenses the .csl list by eliminating CSFs that make no contri-   *
*   bution above the threshold THRESH.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
*                        LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*   Updated by Anders Ynnerman            Last revision: 31 Jan 1994   *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

Cww      INTEGER PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(1))
      CHARACTER*256 FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
      CHARACTER*2 NH
*
      EXTERNAL IQ,ISPAR,ITJPO,JCUP,JQS
*
      POINTER (PINDEX,INDEX(1))
*
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PIATJP,IATJPO(1)),(PIASPA,IASPAR(1))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
      CHARACTER*300 LINE1, LINE2, LINE3
      CHARACTER*500 string(5)
      CHARACTER*24 NAME
*
*   Determine the cutoff criterion
*
      PRINT *, 'What is the value below which an'
      PRINT *, ' eigenvector component is to be'
      PRINT *, ' neglected?'
      READ *, THRESH
*
*    Establish the list of CSFs
*
      CALL ALLOC (PINDEX,NCF,4)
      NSHORT = 0
      DO 2 I = 1,NCF
         DO 1 J = 1,NVEC
            IF (ABS (EVEC(I+(J-1)*NCF)) .GT. THRESH) THEN
               NSHORT = NSHORT+1
               INDEX(NSHORT) = I
               GOTO 2
            ENDIF
    1    CONTINUE
    2 CONTINUE
*
*   Open the output  .csl file; it is FORMATTED; it must not exist
*
      OPEN (UNIT = 21,FILE='rcsl.out',FORM='FORMATTED',
     :     STATUS='UNKNOWN')
      OPEN (UNIT = 19,FILE=TRIM(NAME)//'.c',FORM='FORMATTED',
     :     STATUS='UNKNOWN')
*
*   Determine the common core subshells; write out the list;
*   determine the peel subshells; write out the list; these
*   data form the first part of the header of the .csl file;
*   one additional line forms the remainder of the header of
*   the  .csl file
*
c     NCORE = 0
c     DO 5 J = 1,NW
c        IQFULL = NKJ(J)+1
c        DO 4 I = 1,NSHORT
c           IF (IQ (J,INDEX(I)) .NE. IQFULL) GOTO 6
c   4    CONTINUE
c        NCORE = NCORE+1
c   5 CONTINUE
c     stop
*
*
*   Now write out all CSFs in the reference list
*
    6 DO I=1,5
      READ (19, '(A)') STRING(I)
      WRITE (21, '(A)') TRIM(STRING(I))
      ENDDO
      NCF=0
      DO 7 I = 1,NSHORT
         ICSF = INDEX(I)
   10 READ (19, '(A)') LINE1
      READ (19, '(A)') LINE2
      READ (19, '(A)') LINE3
      NCF=NCF+1
      if(ICSF.EQ.NCF) THEN 
      WRITE(21,'(A)') TRIM(LINE1)
      WRITE(21,'(A)') TRIM(LINE2)
      WRITE(21,'(A)') TRIM(LINE3)
        GO TO 7
      ELSE
      GOTO 10
      ENDIF
c        CALL OUTCSF (ICSF,NCORE,NW,IQ,ISPAR,ITJPO,JCUP,JQS)
    7 CONTINUE
*
      CALL CONVRT (NW,FILNAM,LENTH)
      PRINT *, FILNAM(1:LENTH)//' relativistic subshells;'
      CALL CONVRT (NSHORT,FILNAM,LENTH)
      PRINT *, FILNAM(1:LENTH)//' relativistic CSFs.'
*
      CLOSE (21)
      CLOSE (19)
*
      RETURN
*
  300 FORMAT (120(1X,1I2,1A2))
*
      END
