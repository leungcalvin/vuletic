************************************************************************
*                                                                      *
      SUBROUTINE SETDBG
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 10 Dec 1992   *
*                                                                      *
************************************************************************
!$Id: setdbg.f,v 1.1 2003/09/30 05:51:40 georgio Exp $
!$Log: setdbg.f,v $
!Revision 1.1  2003/09/30 05:51:40  georgio
!
!added
!
!Revision 1.2  1997/03/04 00:37:08  xhh
!I/O stuff.
!
! Short output lines joined
! PRINT *, --> WRITE(istde,*)
! Unit number 99 checked before use
! XHH 1997.01.21

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
      LOGICAL open99
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
         DEFNAM = 'rscf92.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         WRITE(istde,*) 'File  rscf92.dbg  will be created as the'
     &                , ' RSCF92 DeBuG Printout File;'
         WRITE(istde,*) 'enter another file name if this is not'
     &                , ' acceptable; null otherwise:'
         READ (*,'(A)') FILNAM
*
         IF ( LEN_TRIM(FILNAM) .EQ. 0) FILNAM = DEFNAM

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! It is not easy to change this 99 to a variable since it is to be used
! elsewhere outside this subroutine. One thing can be done for safety
! is to add an inquire statement and issue a warning before stoppng
! the program if the unit had been openned before.
! XHH 1997.01.21
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         INQUIRE (UNIT=99, OPENED=open99)
         IF (open99) THEN
            WRITE(istde,*) 'Unit 99 was opened from somewhere'
            STOP
         ENDIF
         
    4    CALL OPENFL (99,FILNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
    5       WRITE(istde,*) 'Enter a name for the RSCF92 DeBuG Printout'
     &, ' file that is to be created:'
            READ (*,'(A)') FILNAM
            IF ( LEN_TRIM(FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         WRITE(istde,*) 'Print out the machine constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(1) = .TRUE.
         WRITE(istde,*) 'Print out the physical constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(2) = .TRUE.
*
         WRITE(istde,*) 'Printout from FNDBLK?'
         YES = GETYN ()
         IF (YES) LDBPG(3) = .TRUE.
         WRITE(istde,*) 'Print out the Hamiltonian matrix?'
         YES = GETYN ()
         IF (YES) LDBPG(4) = .TRUE.
         WRITE(istde,*) 'Print out the eigenvectors?'
         YES = GETYN ()
         IF (YES) LDBPG(5) = .TRUE.
*
*   Set options for printout from radial modules
*
         WRITE(istde,*) 'Printout from RADGRD?'
         YES = GETYN ()
         IF (YES) LDBPR(1) = .TRUE.
         WRITE(istde,*) 'Printout from NUCPOT?'
         YES = GETYN ()
         IF (YES) LDBPR(2) = .TRUE.
         WRITE(istde,*) 'Printout from LODRWF?'
         YES = GETYN ()
         IF (YES) LDBPR(3) = .TRUE.
         WRITE(istde,*) 'Print out I(ab) integrals?'
         YES = GETYN ()
         IF (YES) LDBPR(4) = .TRUE.
         WRITE(istde,*) 'Print out Slater integrals?'
         YES = GETYN ()
         IF (YES) LDBPR(10) = .TRUE.
         WRITE(istde,*) 'Make summary printout on progress'
     &, ' of each iteration in SOLVE?'
         YES = GETYN ()
         IF (YES) LDBPR(22) = .TRUE.
         WRITE(istde,*) 'Tabulate and make printer plots'
     &, ' of subshell radial functions on'
     &, ' each iteration in SOLVE?'
         YES = GETYN ()
         IF (YES) LDBPR(23) = .TRUE.
         WRITE(istde,*) 'Tabulate and make printer plots'
     &, ' of subshell radial functions'
     &, ' after each SCF cycle?'
         YES = GETYN ()
         IF (YES) LDBPR(24) = .TRUE.
         WRITE(istde,*) 'Tabulate and make printer plots'
     &, ' of subshell radial functions on'
     &, ' convergence?'
         YES = GETYN ()
         IF (YES) LDBPR(25) = .TRUE.
         WRITE(istde,*) 'List compositions of exchange'
     &, ' potentials?'
         YES = GETYN ()
         IF (YES) LDBPR(27) = .TRUE.
         WRITE(istde,*) 'Tabulate and make printer plots'
     &, ' of exchange potentials?'
         YES = GETYN ()
         IF (YES) LDBPR(28) = .TRUE.
         WRITE(istde,*) 'List compositions of direct'
     &, ' potentials?'
         YES = GETYN ()
         IF (YES) LDBPR(29) = .TRUE.
         WRITE(istde,*) 'Tabulate and make printer plots'
     &, ' of direct potentials?'
         YES = GETYN ()
         IF (YES) LDBPR(30) = .TRUE.
*
*   Set options for printout of angular coefficients
*
         WRITE(istde,*) ' Printout from LODCSL?'
         YES = GETYN ()
         IF (YES) LDBPA(1) = .TRUE.
         WRITE(istde,*) ' Print out T coefficients?'
         YES = GETYN ()
         IF (YES) LDBPA(2) = .TRUE.
         WRITE(istde,*) ' Print out V coefficients?'
         YES = GETYN ()
         IF (YES) LDBPA(3) = .TRUE.
*
      ENDIF
*
      RETURN
      END
