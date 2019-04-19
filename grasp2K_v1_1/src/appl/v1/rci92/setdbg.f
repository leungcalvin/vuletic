************************************************************************
*                                                                      *
      SUBROUTINE SETDBG
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, LENGTH, OPENFL.                        *
*                                                                      *
*   Written by Farid A Parpia               Last update: 21 Oct 1992   *
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
      DO 3 I = 1,30
         LDBPR(I) = .FALSE.
    3 CONTINUE
*
      IF (NDEF.EQ.0) THEN
         RETURN
      ENDIF

      WRITE(istde,*) 'Generate debug printout?'
      YES = GETYN ()
      IF (YES) THEN
*
*   The  .dbg  file is formatted; open it on unit 99
*
         DEFNAM = 'rci92.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         WRITE(istde,*) 'Enter name for the printout (null==rci92.dbg)'
         READ (*,'(A)') FILNAM
*
         IF ( LEN_TRIM(FILNAM) .EQ. 0) FILNAM = DEFNAM
*
    4    CALL OPENFL (99,FILNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
    5       WRITE(istde,*) 'Input wrong, re-enter the name'
            READ (*,'(A)') FILNAM
            IF ( LEN_TRIM(FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         WRITE(istde,*) ' Print out the machine constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(1) = .TRUE.
         WRITE(istde,*) ' Print out the physical constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(2) = .TRUE.
*
         WRITE(istde,*) ' Printout from FNDBLK?'
         YES = GETYN ()
         IF (YES) LDBPG(3) = .TRUE.
         WRITE(istde,*) ' Print out the Hamiltonian matrix?'
         YES = GETYN ()
         IF (YES) LDBPG(4) = .TRUE.
*
*   Set options for radial modules
*
         WRITE(istde,*) ' Printout from radial modules?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) ' Printout from RADGRD?'
            YES = GETYN ()
            IF (YES) LDBPR(1) = .TRUE.
            WRITE(istde,*) ' Printout from NUCPOT?'
            YES = GETYN ()
            IF (YES) LDBPR(2) = .TRUE.
            WRITE(istde,*) ' Printout from LODRWF?'
            YES = GETYN ()
            IF (YES) LDBPR(3) = .TRUE.
*
            WRITE(istde,*) ' Print out one-electron integrals?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE(istde,*) ' Print out I(ab) integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(4) = .TRUE.
               WRITE(istde,*) ' Print out kinetic energy integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(5) = .TRUE.
               WRITE(istde,*) ' Print out Vinti integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(6) = .TRUE.
               WRITE(istde,*) ' Printout from BESSEL?'
               YES = GETYN ()
               IF (YES) LDBPR(7) = .TRUE.
               WRITE(istde,*) ' Printout from VAC4?'
               YES = GETYN ()
               IF (YES) LDBPR(8) = .TRUE.
               WRITE(istde,*)' Print out vacuum polarisation integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(9) = .TRUE.
            ENDIF
*
            WRITE(istde,*) ' Print out two-electron integrals?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE(istde,*) ' Print out Slater integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(10) = .TRUE.
               WRITE(istde,*) ' Print out Breit integrals?'
               YES = GETYN ()
               IF (YES) LDBPR(11) = .TRUE.
            ENDIF
         ENDIF
*
*   Set options for angular modules
*
         WRITE(istde,*) ' Printout from angular modules?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) ' Printout from LODCSL?'
            YES = GETYN ()
            IF (YES) LDBPA(1) = .TRUE.
            WRITE(istde,*) ' Print out T coefficients?'
            YES = GETYN ()
            IF (YES) LDBPA(2) = .TRUE.
            WRITE(istde,*) ' Print out Coulomb V coefficients?'
            YES = GETYN ()
            IF (YES) LDBPA(3) = .TRUE.
            WRITE(istde,*) ' Print out Breit V coefficients?'
            YES = GETYN ()
            IF (YES) LDBPA(4) = .TRUE.
         ENDIF
*
      ENDIF
*
      RETURN
      END
