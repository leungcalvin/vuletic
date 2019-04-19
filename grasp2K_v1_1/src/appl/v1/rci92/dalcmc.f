************************************************************************
*                                                                      *
      SUBROUTINE DALCMC
*                                                                      *
*   This  subprogram deallocates storage for certain arrays that are   *
*   used by the MCDF section of GRASP2.                                *
*                                                                      *
*   Subprogram(s) called: DALLOC.                                      *
*                                                                      *
*   Written by Farid A. Parpia.           Last revision: 05 Aug 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
Cww      INTEGER PNEVAL,
Cww     :        PNEVEC,
Cww     :        PNIVEC,
Cww     :        PIATJP,PIASPA
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PNEVEC,EVECDUMMY)
      POINTER (PNIVEC,IVECDUMMY)
      POINTER (PIATJP,IATJPDUMMY)
      POINTER (PIASPA,IASPADUMMY)
*
      COMMON/EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
*
*   Deallocate storage for arrays
*
      CALL DALLOC (PNEVAL)
      CALL DALLOC (PNEVEC)
*
      CALL DALLOC (PNIVEC)
*
      CALL DALLOC (PIATJP)
      CALL DALLOC (PIASPA)
*
      RETURN
      END
