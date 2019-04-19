************************************************************************
*                                                                      *
      SUBROUTINE GETINF
*                                                                      *
*   Interactively determines data governing the generation of MCP co-  *
*   efficients.                                                        *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, GETYN.                                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG,GETYN,LFORDR,YES
      CHARACTER*20 CNUM
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEFAULT/NDEF
     :      /FOPARM/ICCUT
     :      /MCPA/KMAX
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
!XHH I/O units
      COMMON/iounit/istdi,istdo,istde
*
*   Determine the physical effects specifications
*
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'Generate MCP coefficients only for'
     & , ' diagonal matrix elements? '
         WRITE(istde,*) '(This is appropriate to an (E)AL calculation):'
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF 
      IF (YES) THEN
         DIAG = .TRUE.
         LFORDR = .FALSE.
         ICCUT = 0
      ELSE
         DIAG = .FALSE.
         IF (NDEF.NE.0) THEN
            WRITE(istde,*) 'Treat contributions of some CSFs'
     &,              ' as first-order perturbations?'
            YES = GETYN ()
         ELSE
            YES = .FALSE.
         ENDIF
         IF (YES) THEN
            LFORDR = .TRUE.
            WRITE(istde,*) 'The contribution of CSFs 1 -- ICCUT will'
     &,              ' be treated variationally;'
            WRITE(istde,*) 'the remainder perturbatively; enter ICCUT:'
    1       READ *, ICCUT
            IF ((ICCUT .LE. 1) .OR. (ICCUT .GE. NCF)) THEN
               CALL CONVRT (NCF,CNUM,LENTH)
               WRITE(istde,*) 'GETINF: ICCUT must be greater than 1'
     &,                 ' and less than '//CNUM(1:LENTH)//';'
               WRITE(istde,*) ' please reenter ICCUT:'
               GOTO 1
            ENDIF
         ELSE
            LFORDR = .FALSE.
            ICCUT = 0
         ENDIF
      ENDIF
*
*   Add the second and third records to the file headers; the
*   first is written in SETRES
*
      DO 2 K = 30,32+KMAX
         WRITE (K) NELEC,NCF,NW
         WRITE (K) DIAG,ICCUT,LFORDR
    2 CONTINUE
*
      RETURN
      END
