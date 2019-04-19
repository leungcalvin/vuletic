************************************************************************
*                                                                      *
      SUBROUTINE DNICMV (N,M,B,C)
*                                                                      *
*   Matrix-matrix product: C = AB.  The lower triangle of the  (NxN)   *
*   matrix is assumed available in packed form in the array EMT. The   *
*   matrices B and C are (NxM).                                        *
*                                                                      *
*   This is an adaptation of  Andreas Stathopulos'  routine  TRSBMV,   *
*   and is specific to GRASP2 derivatives.                             *
*                                                                      *
*   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
*               [BLAS]: DAXPY/SAXPY, DDOT/SDOT.                        *
*                                                                      *
*   F A Parpia and A Stathopoulos         Last revision: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
Cww      INTEGER PIENDC,PNELBL,PNEVAL,PNIROW
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNELBL,ELBLDUMMY)
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PNIROW,IROWDUMMY)
*
      POINTER (PIBLOC,IBLOCK(1))
      POINTER (PNTEMT,EMT(1))
*
      COMMON/EIGVAL/EAV,PNEVAL
     :      /HBLOCK/NBLOCK,PIBLOC,PNELBL
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /WCHBLK/JBLOCK
*
      DIMENSION B(N,M),C(N,M)
*
*   Initialise the result matrix; note that this is specific to the
*   data structure of DVDSON --- there is no overdimensioning
*
      CALL DINIT (N*M,0.0D 00,C,1)
*
      ICUR = 1
      DO 2 ICOL = 1,N
         IF (IBLOCK(ICOL) .EQ. JBLOCK) THEN
	    DO 1 IV = 1,M
               C(ICOL,IV) =  C(ICOL,IV)
     :                      +EMT(ICUR)*B(ICOL,IV)
     :              +DDOT (N-ICOL,EMT(ICUR+1),1,
     :                                   B(ICOL+1,IV),1)
               CALL DAXPY (N-ICOL,B(ICOL,IV),EMT(ICUR+1),1,
     :				                    C(ICOL+1,IV),1)
    1       CONTINUE
         ENDIF
         ICUR=ICUR+N-ICOL+1
    2 CONTINUE
*
      RETURN
      END
