************************************************************************
*                                                                      *
      SUBROUTINE SETSDA (NNONZ,LPRINT)
*                                                                      *
*   This routine examines lists                                        *
*                                                                      *
*                               (IC,IR,INDEX)                          *
*                                                                      *
*   to set up the array  IENDC  required by the Davidson eigensolver   *
*   of Stathopoulos and Fischer.                                       *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL DIAG,LFORDR,LPRINT
      CHARACTER*20 CNUM
      CHARACTER*8 SRTLAB
      CHARACTER*3 MCPLAB
*
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
cff  
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/ORB2/NCF,NW,PNTRIQ
*
      REWIND (30)
*
*   Has the file already been analysed? Skip the analysis if so.
*
      READ (30) MCPLAB,SRTLAB
      IF (SRTLAB .EQ. '  SORTED') GOTO 2
*
*   Message
*
      PRINT *, 'Analysing sparse matrix array definition file ...'
*
*   Position the file
*
      READ (30) NELEC,NCF,NW
      READ (30) DIAG,ICCUT,LFORDR
*
*   Allocate storage for IENDC(0:NCF)
*
      CALL ALLOC (PIENDC,NCF+1,4)
      CALL ALLOC (PNIROW,NNONZ,4)
*
*   Analyse data on file 30; set up IENDC and IROW
*
      IEND = 0
      ICLAST = 0
*
      DO 1 I = 1,NNONZ
         READ (30) IC,IROW(I),INDEX
         IF (IC .NE. ICLAST) THEN
            IENDC(ICLAST) = IEND
            ICLAST = IC
         ENDIF
         IEND = INDEX
    1 CONTINUE
      IENDC(NCF) = IEND
*
*   Sorting complete; rewrite the file header
*
      REWIND (30)
      WRITE (30) 'MCP','  SORTED'
      WRITE (30) NELEC,NCF,NW
      WRITE (30) DIAG,ICCUT,LFORDR
*
*   Write the number of elements to file 30
*
      WRITE (30) NNONZ
*
*   Write arrays IENDC and IROW to file 30
*
      WRITE (30) (IENDC(I),I = 0,NCF),(IROW(I),I = 1,NNONZ)
*
*   Deallocate storage
*
      CALL DALLOC (PIENDC)
      CALL DALLOC (PNIROW)
*
*   Completion message
*
      CALL CONVRT (NNONZ,CNUM,LCNUM)
      PRINT *, ' ... analysis complete; '//CNUM(1:LCNUM)//' nonzero'
      NMAX = (NCF*(NCF+1))/2
      CALL CONVRT (NMAX,CNUM,LCNUM)
      PRINT *, ' elements in H(DC); maximum possible: '
     :   //CNUM(1:LCNUM)//';'
*
*   Debug printout
*
    2 IF (LPRINT) THEN
         WRITE (99,300)
         REWIND (30)
         READ (30)
         READ (30)
         READ (30)
         READ (30) NELMNT
         CALL ALLOC (PIENDC,NCF+1,4)
         CALL ALLOC (PNIROW,NELMNT,4)
         WRITE (99,301) NELMNT
         READ (30) (IENDC(I),I = 0,NCF),(IROW(I),I = 1,NELMNT)
         DO 4 IC = 1,NCF
            DO 3 I = IENDC(IC-1)+1,IENDC(IC)
               WRITE (99,302) IC,IROW(I),I
    3       CONTINUE
    4    CONTINUE
         CALL DALLOC (PIENDC)
         CALL DALLOC (PNIROW)
      ENDIF
*
*   Close file 30; it is no longer required by GENMCP
*
      CLOSE (30)
*
      RETURN
*
  300 FORMAT (/'From SETSDA:')
  301 FORMAT (' Number of nonzero elements in H(DC): ',1I4)
  302 FORMAT (' Column ',1I2,', row ',1I2,', sparse matrix index ',1I4)
*
      END
