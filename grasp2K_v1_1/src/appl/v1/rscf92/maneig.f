************************************************************************
*                                                                      *
      SUBROUTINE MANEIG (LPRINT)
! new version developed by Misha
! the old version is renamed temporarily as maneig_old.f
! 97.02.12
*                                                                      *
*   This module  manages the  operation of the  eigensolvers and the   *
*   storage of the eigenpairs.  There are two principal branches:      *
*                                                                      *
*      (1) Matrix of order 1: the trivial case                         *
*      (2) Matrix of order greater than 1: eigenpairs are found        *
*          using DVDSON; this involves up to three steps:              *
*                    (a) The matrix is analysed to determine its       *
*                        block structure (only irreducibe matrices     *
*                        are correctly treated by DVDSON)              *
*                    (b) Eigenpairs are extracted for each block       *
*                    (c) The appropriate eigenpairs are selected and   *
*                        stored                                        *
*                                                                      *
*   We  assume that  the sparse representation  of the matrix  is in   *
*   core.                                                              *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, RALLOC.          *
*               [RSCF92]: FNDBLK, POSNFL, SPICMV.                      *
*               [DVDSON]: DVDSON                                       *
*               [AUXBLAS]: DINIT/SINIT                                 *
*               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL, DSWAP/SSWAP          *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL FIRST,HIEND,LPRINT
*
      EXTERNAL SPICMV
*
      POINTER (PNWORK,WORK(1))
      POINTER (PIWORK,IWORK(1))
      POINTER (PNDIAG,DIAG(1))
      POINTER (PJWORK,JWORK(1))
*
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PIBLOC,IBLOCK(1))
      POINTER (PNELBL,NELBLK(1))
      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
*
      COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /HBLOCK/NBLOCK,PIBLOC,PNELBL
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /SYMA/PIATJP,PIASPA
     :      /WCHBLK/JBLOCK
     :      /WHERE/IMCDF,NREC
*
      DATA FIRST /.TRUE./
*
*   Entry message
*
      WRITE (*,300)
*
      IF (NCF .EQ. 1) THEN
*
*   Matrix of order 1: the trivial case; we assume that the value
*   of EAV is available
*
         CALL ALLOC (PNEVAL,1,8)
         CALL ALLOC (PNEVEC,1,8)
         EVAL(1) = 0.0D 00
         EVEC(1) = 1.0D 00
*
      ELSE
*
*   Use Davidson eigensolver
*
*   Analyse the Hamiltonian matrix to determine its block structure
*
         CALL FNDBLK (LPRINT)
*
*   Store the diagonals in a separate array
*
         CALL ALLOC (PNDIAG,NCF,8)
*
         DO 1 IC = 1,NCF
            DIAG(IC) = EMT(IENDC(IC-1)+1)
    1    CONTINUE
*
*   Allocate storage for workspace; see the header of DVDSON for
*   the expression below; the value of LIM can be reduced to NVECT
*   plus a smaller number if storage is severely constrained
*
         NVECT = NCMAX
         LIM = MIN (NCF,NVECT+20)
c        LWORK = 2*NCF*LIM+LIM*LIM+(NVECT+10)*LIM+NVECT
         LWORK = 2*NCF*LIM+LIM*LIM*2+11*LIM+NVECT
         CALL ALLOC (PNWORK,LWORK,8)
         LIWORK = 6*LIM+NVECT
         CALL ALLOC (PIWORK,LIWORK,4)
*
         MBLOCK = 1
         
*changed by Misha 02/03/97
         CRITE = 1.0D-17
         CRITC = 1.0D-08
         CRITR = 1.0D-08
         ORTHO = max(1D-8,CRITR)
* end of changes

Cww         MAXITR = MAX (NVECT*100,NCF/10)
         MAXITR = MIN (NVECT*100,NCF)
         CALL ALLOC (PJWORK,LIM,4)
*
*   Matrix is reducible; find MIN (NVECT,NELBLK(I)) lowest
*   eigenpairs in each block I
*
         CALL ALLOC (PNEVAL,    NVECT,8)
         CALL ALLOC (PNEVEC,NCF*NVECT,8)
*
Cww  Error
Cww      DMUNGO = 10.0D 00**(TENMAX-1.0D 00)
         DMUNGO = 10.0D 00**(308)
         CALL DINIT (NVECT,DMUNGO,EVAL,1)
*
*         CALL ALLOC (PJWORK,LIM,4)
*
*   Compute the eigenpairs in each block
*
         DO 12 JBLOCK = 1,NBLOCK
*
            NVEX = MIN (NVECT,NELBLK(JBLOCK),14)
            NEND = NCF*NVEX
*
            ILOW = 1
            IHIGH = NVEX
            NIV = NVEX
*
*   Initial estimates for eigenvectors
*
*changed by Misha 02/03/97
	    call INIEST(NCF,NELBLK(JBLOCK),NIV,EMT,IENDC(0),IROW,WORK,
     :                 IBLOCK, JBLOCK)

* iniest looks for eigenvectors of 1000*1000 matrix so there is no need
* to call dvdson if block size <= 1000

	    if (NELBLK(JBLOCK) .GT. 1000 ) then
!misha back	    if (NELBLK(JBLOCK) .GT. 1 ) then
	    write (*,*) 'CALLING DVDSON!!!', maxitr

*
*   Call Davidson eigensolver
*
      print *
      print *
      print *, ' Inside maneig'
      print *
      print *, ' NCF = ', NCF
      print *, ' LIM = ', LIM
      print *, ' ILOW = ', ILOW
      print *, ' IHIGH = ', IHIGH
      print *, ' NIV = ', NIV
      print *, ' MBLOCK = ', MBLOCK
      print *, ' ORTHO = ', ORTHO
      print *, ' MAXITR = ', MAXITR
      print *, ' LWORK = ', LWORK
      print *, ' LIWORK = ', LIWORK
      print *, ' WORK(1) = ', WORK(1)
      print *
      print *
              CALL GDVD (SPICMV,NCF,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)
      print *, ' HIEND = ', HIEND
      print *, ' NLOOPS = ', NLOOPS
      print *, ' NMV = ', NMV
      print *, ' IERR = ', IERR
      print *, ' WORK(1) = ', WORK(1)
      print *
      print *
*
              WRITE (*,301) NLOOPS,NMV
*
              IF (IERR .NE. 0) THEN
                 WRITE (*,302) IERR
*misha                 STOP
              ENDIF
      
           endif
           
*
*   Put the eigenpairs in order, overwriting as necessary
*
            IF (JBLOCK .EQ. 1) THEN
               CALL DCOPY (    NVEX,WORK(NEND+1),1,EVAL,1)
               CALL DCOPY (NCF*NVEX,WORK(     1),1,EVEC,1)
            ELSE
               DO 11 I = 1,NVEX
                  EVALI = WORK(NEND+I)
                  DO 10 J = I,NVECT
                     IF (EVALI .LT. EVAL(J)) THEN
                        DO 9 K = NVECT,J+1,-1
                           EVAL(K) = EVAL(K-1)
                           CALL DCOPY (NCF,
     :                            EVEC((K-2)*NCF+1),1,
     :                            EVEC((K-1)*NCF+1),1)
    9                   CONTINUE
                        EVAL(J) = EVALI
                        CALL DCOPY (NCF,
     :                            WORK((I-1)*NCF+1),1,
     :                            EVEC((J-1)*NCF+1),1)
                        GOTO 11
                     ENDIF
   10             CONTINUE
   11          CONTINUE
            ENDIF
*
   12    CONTINUE
*
*   Deallocate storage
*
         CALL DALLOC (PNTEMT)
         CALL DALLOC (PNDIAG)
         CALL DALLOC (PNWORK)
         CALL DALLOC (PIWORK)
         CALL DALLOC (PJWORK)
         CALL DALLOC (PIBLOC)
         CALL DALLOC (PNELBL)
*
*   Rearrange and reallocate storage for the eigenpairs
*   as necessary
*
         CALL ALLOC (PIWORK,NVECT,4)
         DO 13 I = 1,NVECT
            IWORK(I) = I
   13    CONTINUE
         NVEX = NCMIN
         DO 14 I = 1,NVEX
            IOFSET = ICCMIN(I)
            LOC = IWORK(I)
            IF (IOFSET .NE. LOC) THEN
               CALL DSWAP (1,EVAL(IOFSET),1,
     :                               EVAL(I     ),1)
               IWORK(I) = IWORK(IOFSET)
               IWORK(IOFSET) = LOC
               IOFSET = NCF*(IOFSET-1)
               LOC = NCF*(I-1)
               CALL DSWAP (NCF,EVEC(IOFSET+1),1,
     :                                 EVEC(LOC   +1),1)
            ENDIF
   14    CONTINUE
         CALL DALLOC (PIWORK)
         CALL RALLOC (PNEVAL,    NVECT,    NVEX,8)
         CALL RALLOC (PNEVEC,NCF*NVECT,NCF*NVEX,8)
*
      ENDIF
*
*   Clean up eigenvectors; determine their J/P values
*
      NVEX = NCMIN
      DO 17 J = 1,NVEX
*
*   Find the dominant component of each eigenvector
*
         IOFSET = (J-1)*NCF
*
         AMAX = 0.0D 00
         DO 15 I = 1,NCF
            WA = ABS (EVEC(I+IOFSET))
            IF (WA .GT. AMAX) THEN
               AMAX = WA
               IA = I
            ENDIF
   15    CONTINUE
*
*   Find the angular momentum and parity of the dominant component
*
         ITJPIA = ITJPO (IA)
         ISPAIA = ISPAR (IA)
*
*   For (E)OL calculations, write out a warning message if these
*   have changed since the last iteration
*
         IF (.NOT. FIRST) THEN
            IF ( (ITJPIA .NE. IATJPO(J)) .OR.
     :           (ISPAIA .NE. IASPAR(J)) ) THEN
               WRITE (*,303) J,IATJPO(J),IASPAR(J),
     :                         ITJPIA   ,ISPAIA
            ENDIF
         ENDIF
*
*   Assign these to the ASF
*
         IATJPO(J) = ITJPIA
         IASPAR(J) = ISPAIA
*
*   Remove any spurious component from the eigenvectors; compute
*   the normalization factor at the same time
*
         SUM = 0.0D 00
         DO 16 I = 1,NCF
            IF ((ITJPO (I) .NE. ITJPIA) .OR.
     :          (ISPAR (I) .NE. ISPAIA)) THEN
               EVEC(I+IOFSET) = 0.0D 00
            ELSE
               SUM = SUM+EVEC(I+IOFSET)**2
            ENDIF
   16    CONTINUE
         DNFAC = 1.0D 00/SQRT (SUM)
*
*   Redefine eigenvectors so that the dominant component
*   is positive
*
         IF (EVEC(IA+IOFSET) .LT. 0.0D 00)
     :      DNFAC = -DNFAC
*
*   Perform replacement with normalization and/or inversion
*
         CALL DSCAL (NCF,DNFAC,EVEC(IOFSET+1),1)
*
   17 CONTINUE
*
      FIRST = .FALSE.
*
      RETURN
*
  300 FORMAT (/'MANEIG ...'/)
  301 FORMAT ('DVDSON: ',1I3,' loops; ',
     :                   1I3,' matrix-vector multiplies.')
  302 FORMAT (' Returned from DVDSON with IERR = ',1I4)
  303 FORMAT (/' ***** WARNING *****'
     :       //' The angular momentum and parity of level ',1I2,
     :         ' have changed:'
     :        /' Last iteration: (2J+1) = ',1I2,', parity = ',1I2,';'
     :        /' this iteration: (2J+1) = ',1I2,', parity = ',1I2,'.')
*
      END
