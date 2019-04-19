************************************************************************
*                                                                      *
      SUBROUTINE SETRES(NAME)
*                                                                      *
*   Open, check, load data from the  .res  file.                       *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, LENGTH, OPENFL.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL FOUND,GETYN,RESTRT,YES
      CHARACTER*256 FILNAM
      CHARACTER*24 NAME
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*6 R92RES
      CHARACTER*3 STATUS
      COMMON/DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde
*
*   File  rci.res  is UNFORMATTED
*
      DEFNAM = 'rci.res'
      FORM = 'UNFORMATTED'
*
*   Determine if this is a restart
*
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'Restarting RCI92 ?'
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF
*
      IF (YES) THEN
*
*   Assume that the computation can be restarted
*
         RESTRT = .TRUE.
         STATUS = 'OLD'
*
*   Look for `rci.res'
*
         INQUIRE (FILE = DEFNAM, EXIST = FOUND)
*
         IF (FOUND) THEN
*
             FILNAM = DEFNAM
*
         ELSE
*
    1       WRITE(istde,*) 'Error when opening rci.res'
            STOP
*
         ENDIF
*
      ELSE
*
*   Not a restart
*
         RESTRT = .FALSE.
         STATUS = 'NEW'
*
         FILNAM = DEFNAM
*
      ENDIF
*
      CALL OPENFL (26,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         IF (RESTRT) THEN
            STOP
         ELSE
            STOP
         ENDIF
      ENDIF
*
      IF (RESTRT) THEN
*
*   Check the file if restarting
*
         READ (26,IOSTAT = IOS) R92RES
         IF ((IOS .NE. 0) .OR. (R92RES .NE. 'R92RES')) THEN
            WRITE(istde,*) 'Not an RCI92 REStart File;'
            CLOSE (26)
            STOP
         ENDIF
*
*   Read and check restart information
*
         CALL LODRES
*
      ELSE
*
*   Write the file header
*
         WRITE (26) 'R92RES'
*
*   Generate the first part of the  .res  file
*
         CALL GETCID(NAME)
*
      ENDIF
*
      RETURN
      END
