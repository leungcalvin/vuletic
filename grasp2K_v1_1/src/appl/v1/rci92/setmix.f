************************************************************************
*                                                                      *
      SUBROUTINE SETMIX(NAME)
*                                                                      *
*   Opens the  .mix  file on stream 25; writes a header to this file;  *
*   calls LODMIX to interactively determine the eigenpairs required.   *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*               [RCI92]: LODMIX.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*24 NAME
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
*   File  grasp92.mix  is UNFORMATTED
*
      K = INDEX(NAME,' ')
      FILNAM = NAME(1:K-1)//'.cm'
      FORM = 'UNFORMATTED'
*
      STATUS = 'NEW'
*
      CALL OPENFL (25,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error when opening',FILNAM
         STOP
      ENDIF
*
      WRITE (25) 'G92MIX'
*
      CALL LODMIX
*
      RETURN
      END
