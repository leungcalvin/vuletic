************************************************************************
*                                                                      *
      SUBROUTINE SETRWG
*                                                                      *
*   Open and write a header to the output .rwf  file.                  *
*                                                                      *
*   Call(s) to: [LIB92]: DALLOC, LENGTH, OPENFL.                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
! The same file will be openned again in orbout.f where the radial 
! functions are to be written to the disk. 
! Note if the filename 'rwfn.out' below is changed, then must be
! the filename in orbout.f
! The call to this routine from rscf92.f should better be removed
! Yes, removed.
! XHH 1997.02.14

      CHARACTER*256 FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
      COMMON/iounit/istdi,istdo,istde
*
*   File  .rwf  file is UNFORMATTED; it is to be created
*
      FORM = 'UNFORMATTED'
      STATUS = 'NEW'
      FILNAM = 'rwfn.out'
*
      CALL OPENFL (23,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) THEN
         WRITE(istde,*) 'Error when opening rwfn.out'
         STOP
      ENDIF
*
*   Write the file header
*
      WRITE (23) 'G92RWF'
*
      RETURN
      END
