************************************************************************
*                                                                      *
      SUBROUTINE SETDBG
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
*                                                                      *
*   Written by Farid A Parpia               Last update: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN,LDBPA,LDBPG,LDBPR,YES
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEBUGR/LDBPR(30)
     :      /DEFAULT/NDEF
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
      DO 3 I = 1,30
         LDBPR(I) = .FALSE.
    3 CONTINUE
*
      IF (NDEF.EQ.0) THEN
         RETURN
      ENDIF

      PRINT *, 'Generate debug printout?'
      YES = GETYN ()
      IF (YES) THEN
*
*   The  .dbg  file is formatted; open it on unit 99
*
         DEFNAM = 'oscl92.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         PRINT *, 'File  oscl92.dbg  will be created as the'
         PRINT *, ' OSCL92 DeBuG Printout File; enter another'
         PRINT *, ' file name if this is not acceptable;'
         PRINT *, ' null otherwise:'
         READ (*,'(A)') FILNAM
*
         IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = DEFNAM
*
    4    CALL OPENFL (99,FILNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
    5       PRINT *, 'Enter a name for the OSCL92 DeBuG Printout'
            PRINT *, ' file that is to be created:'
            READ (*,'(A)') FILNAM
            IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         PRINT *, 'Print out the machine constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(1) = .TRUE.
         PRINT *, 'Print out the physical constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(2) = .TRUE.
*
*   Set options for radial modules
*
         PRINT *, 'Printout from radial modules?'
         YES = GETYN ()
         IF (YES) THEN
            PRINT *, 'Printout from RADGRD?'
            YES = GETYN ()
            IF (YES) LDBPR(1) = .TRUE.
            PRINT *, 'Printout from LODRWF?'
            YES = GETYN ()
            IF (YES) LDBPR(3) = .TRUE.
*
            PRINT *, 'Printout from radial modules?'
            YES = GETYN ()
            IF (YES) THEN
*
               PRINT *, 'Print M integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(12) = .TRUE.
               PRINT *, 'Print integrands of M integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(13) = .TRUE.
               PRINT *, 'Print I and J  integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(14) = .TRUE.
               PRINT *, 'Print integrands of I and J integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(15) = .TRUE.
               PRINT *, 'Print Bessel functions?'
               YES = GETYN ()
               IF (YES) LDBPR(16) = .TRUE.
*
            ENDIF
*
         ENDIF
*
*   Set options for angular modules
*
         PRINT *, 'Printout from angular modules?'
         YES = GETYN ()
         IF (YES) THEN
            PRINT *, 'Printout from LODCSL?'
            YES = GETYN ()
            IF (YES) LDBPA(1) = .TRUE.
            PRINT *, 'Print out T coefficients?'
            YES = GETYN ()
            IF (YES) LDBPA(2) = .TRUE.
         ENDIF
*
      ENDIF
*
      RETURN
      END
