      subroutine setsum(NAME)
*                                                                      *
*   Open the  .sum  file on stream 24.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*24 NAME
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      COMMON/iounit/istdi,istdo,istde
*
*   File  rci92.sum  is FORMATTED
*
      K = INDEX(NAME,' ')
      FILNAM = NAME(1:K-1)//'.csum'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
    1 CALL OPENFL (24,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(istde,*) 'Error when opening',FILNAM
         STOP
      ENDIF
*
      RETURN
      END
