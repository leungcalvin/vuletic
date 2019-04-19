************************************************************************
*                                                                      *
      SUBROUTINE SPODMV (N,M,B,C)
*                                                                      *
*   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
*   lower triangle of the  (NxN)  matrix  A  is assumed available in   *
*   COMMON/HMAT/.                                                      *
*                                                                      *
*   This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,   *
*   and is specific to GRASP2 derivatives.                             *
*                                                                      *
*   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
*               [SPBLAS]: DAXPYI/SAXPYI, DDOTI/SDOTI.                  *
*                                                                      *
*   F A Parpia and A Stathopoulos         Last revision: 13 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
Cww      INTEGER PNELBL,PNEVAL,PIENDC
      POINTER (PNELBL,ELBLDUMMY)
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PIENDC,ENDCDUMMY)
*
      POINTER (PIBLOC,IBLOCK(1))
      POINTER (PNTEMT,EMT(1))
      POINTER (PNIROW,IROW(1))
*
      COMMON/EIGVAL/EAV,PNEVAL
      COMMON/HBLOCK/NBLOCK,PIBLOC,PNELBL
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /WCHBLK/JBLOCK
     :      /WHERE/IMCDF,NREC
*
      DIMENSION B(N,M),C(N,M)
*
*   Initialise the result matrix; note that this is specific to the
*   data structure of DVDSON
*
      CALL DINIT (N*M,0.0D 00,C,1)
*
*   Initialise the file on which the matrix is stored
*
      CALL POSNFL (IMCDF,NREC)
*
      DO 2 ICOL = 1,N
         IF (IBLOCK(ICOL) .EQ. JBLOCK) THEN
            READ (IMCDF) NELC,ELSTO,(EMT(IR),IR = 1,NELC),
     :                             (IROW(IR),IR = 1,NELC)
	    DO 1 IV = 1,M
               DIAG =  C(ICOL,IV)
     :                      +(EMT(1)-EAV)*B(ICOL,IV)
               CALL DMERGE (NELC-1,B(1,IV),C(1,IV),
     :                     IROW(2),EMT(2),B(ICOL,IV),DL)
               C(ICOL,IV) = DIAG + DL
    1       CONTINUE
         ELSE
            READ(IMCDF)
         ENDIF
    2 CONTINUE
*
      RETURN
      END
