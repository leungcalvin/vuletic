************************************************************************
*                                                                      *
      SUBROUTINE PRTRSL
*                                                                      *
*   Prints  subshells that  are to be varied in the order defined by   *
*   IORDER.                                                            *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
! PRINT *, --> WRITE(istde,*)
! XHH 1997.01.22

Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LFIX
      CHARACTER*80 RECORD
      CHARACTER*2 CNUM,NH
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/FIXD/NFIX,LFIX(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORBA/IORDER(NNNW)
      COMMON/iounit/istdi,istdo,istde
*
*   Print a list of subshell radial wavefunctions that remain to
*   be estimated; this list is no more than 80 characters wide
*
      IEND = 0
      DO 1 I = 1,NW
         LOC = IORDER(I)
         IF (.NOT. LFIX(LOC)) THEN
            IF (IEND .GT. 75) THEN
               WRITE(istde,*) RECORD(1:IEND)
               IEND = 0
            ENDIF
            IF (IEND .GT. 0) THEN
               IBEG = IEND+1
               IEND = IBEG
               RECORD(IBEG:IEND) = ' '
            ENDIF
            IBEG = IEND+1
            CALL CONVRT (NP(LOC),CNUM,LENTH)
            IEND = IBEG+LENTH-1
            RECORD(IBEG:IEND) = CNUM(1:LENTH)
            IF (NAK(LOC) .LT. 0) THEN
               LENTH = 1
            ELSE
               LENTH = 2
            ENDIF
            IBEG = IEND+1
            IEND = IBEG+LENTH-1
            RECORD(IBEG:IEND) = NH(LOC)(1:LENTH)
         ENDIF
    1 CONTINUE
      IF (IEND .GT. 1) WRITE(istde,*) RECORD(1:IEND)
*
      RETURN
      END
