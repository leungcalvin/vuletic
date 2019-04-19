************************************************************************
*                                                                      *
      SUBROUTINE SETSUM
*                                                                      *
*   Open the  .sum  file on stream 24.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
!$Id: setsum.f,v 1.1 2003/09/30 05:51:40 georgio Exp $
!$Log: setsum.f,v $
!Revision 1.1  2003/09/30 05:51:40  georgio
!
!added
!
!Revision 1.2  1997/03/04 00:43:05  xhh
!I/O style.
!
! Safety open unit 24 added
! PRINT *, --> WRITE(istde,*)
! XHH 1997.01.21

      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      COMMON/iounit/istdi,istdo,istde
      LOGICAL open24

*
*   File  rscf92.sum  is FORMATTED
*
      DEFNAM = 'rscf.sum'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      FILNAM = DEFNAM
*

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! It is not easy to change this 24 to a variable since it is to be used
! elsewhere outside this subroutine. One thing can be done for safety
! is to add an inquire statement and issue a warning before stoppng
! the program if the unit had been openned before.
! XHH 1997.01.21
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      INQUIRE (UNIT=24, OPENED=open24)
      IF (open24) THEN
         WRITE(istde,*) 'Unit 24 was opened from somewhere'
         STOP
      ENDIF
         
      CALL OPENFL (24,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(istde,*) 'Error when opening rscf.sum'
         STOP
      ENDIF
*
      RETURN
      END
