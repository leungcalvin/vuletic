************************************************************************
*                                                                      *
      SUBROUTINE SETDBG (DEBUG)
*                                                                      *
*   This routine opens the  .dbg  file and sets the arrays that con-   *
*   trol debug printout from the GENMCP program.                       *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
*                                                                      *
*   Written by Farid A Parpia               Last update: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL DEBUG,GETYN,LDBPA,LDBPG,YES
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEFAULT/NDEF

!XHH I/O units
      COMMON/iounit/istdi,istdo,istde
*
*   Initialise the arrays that control the debug printout
*
      DO 1 I = 1,5
         LDBPA(I) = .FALSE.
    1 CONTINUE
*
      DO 2 I = 1,5
         LDBPG(I) = .FALSE.
    2 CONTINUE
*
      IF (NDEF.EQ.0) THEN
         RETURN
      ENDIF

      WRITE(istde,*) 'Generate debug printout?'
      YES = GETYN ()
      IF (YES) THEN
*
         DEBUG = .TRUE.
*
*   The  .dbg  file is formatted; open it on unit 99
*
         DEFNAM = 'genmcp.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         WRITE(istde,*) 'File  genmcp.dbg  will be created as the'
     &,                 ' GENMCP DeBuG Printout File; '
         WRITE(istde,*) ' enter another file name if this is not '
     &,                 'acceptable; null otherwise:'
         READ (*,'(A)') FILNAM
*
         IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = DEFNAM
*
    4    CALL OPENFL (99,FILNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
    5       WRITE(istde,*) 'Enter a name for the GENMCP DeBuG Printout'
     &,                    ' file that is to be created:'
            READ (*,'(A)') FILNAM
            IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         WRITE(istde,*) ' Print out the machine constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(1) = .TRUE.
*
*   Set options for angular modules
*
         WRITE(istde,*) ' Printout from LODCSL?'
         YES = GETYN ()
         IF (YES) LDBPA(1) = .TRUE.
         WRITE(istde,*) ' Print out T coefficients?'
         YES = GETYN ()
         IF (YES) LDBPA(2) = .TRUE.
         WRITE(istde,*) ' Print out Coulomb V coefficients?'
         YES = GETYN ()
         IF (YES) LDBPA(3) = .TRUE.
         WRITE(istde,*) ' Print out sparse matrix definition arrays?'
         YES = GETYN ()
         IF (YES) LDBPA(4) = .TRUE.
*
      ELSE
*
         DEBUG = .FALSE.
*
      ENDIF
*
      RETURN
      END
