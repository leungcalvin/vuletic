************************************************************************
*                                                                      *
      SUBROUTINE SETMIX
*                                                                      *
*   Opens the  .mix  file on stream 25; writes a header to this file.  *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /ORB2/NCF,NW,PNTRIQ
*
*   File  grasp92.mix  is UNFORMATTED
*
      DEFNAM = 'rmix.out'
      FORM = 'UNFORMATTED'
*
      STATUS = 'NEW'
*
      FILNAM = DEFNAM
*
      CALL OPENFL (25,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error when opening rmix.out'
         STOP
      ENDIF
*
*   Write the file header
*
      WRITE (25) 'G92MIX'
      WRITE (25) NELEC,NCF,NW
      WRITE (25) NCMIN
      WRITE (25) (ICCMIN(I),I = 1,NCMIN)
*
      RETURN
      END
