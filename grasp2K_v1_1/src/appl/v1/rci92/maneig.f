************************************************************************
*                                                                      *
      SUBROUTINE MANEIG (LPRINT)
*                                                                      *
*   This module  manages the  operation of the  eigensolvers and the   *
*   storage of the eigenpairs.  There are two principal branches:      *
*                                                                      *
*      (1) Matrix of order 1: the trivial case                         *
*      (2) Matrix of order greater than 1: there are two branches      *
*             (i) Matrices of order less than or equal to IOLPCK:      *
*                 eigenpairs are found using LAPACK SUBROUTINEs        *
*            (ii) Matrices of order greater than IOLPCK: eigenpairs    *
*                 are found using DVDSON; this involves up to three    *
*                 steps:                                               *
*                    (a) The matrix is analysed to determine its       *
*                        block structure (only irreducibe matrices     *
*                        are correctly treated by DVDSON)              *
*                    (b) Eigenpairs are extracted for each block       *
*                    (c) The appropriate eigenpairs are selected and   *
*                        stored                                        *
*                 Different methods of storage and different           *
*                 versions of the matrix-vector multiply are used      *
*                 depending upon the order and density of the matrix   *
*   
*   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, POSNFL,          *
*                        RALLOC.                                       *
*               [RCI92]: DNICMV, FNDBLK, SPICMV, SPODMV.               *
*               [DVDSON]: DVDSON.                                      *
*               [AUXBLAS]: DINIT/SINIT.                                *
*               [BLAS]: DCOPY/SCOPY, DSWAP/SSWAP.                      *
*               [LAPACK]: DSPEVX/SSPEVX.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
*   Modified Misha Saparov                                  Feb 1997   *
*            Charlotte F. Fischer                           May 1997   *
*                 Except for the disk version, all matrices have       *
*                 diagonals  shifted by EAV                            *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL HIEND,LDISC,LPRINT,SPARSE
      CHARACTER*8 CNUM
*
      EXTERNAL DNICMV,SPICMV,SPODMV
*
      POINTER (PONTRW,W(1))
      POINTER (PONTRZ,Z(1))
      POINTER (PNWORK,WORK(1))
      POINTER (PIWORK,IWORK(1))
      POINTER (PIFAIL,IFAIL(1))
      POINTER (PNDIAG,DIAG(1))
      POINTER (PJWORK,JWORK(1))
*
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PIBLOC,IBLOCK(1))
      POINTER (PNELBL,NELBLK(1))
      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
*
      COMMON/EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /HBLOCK/NBLOCK,PIBLOC,PNELBL
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /WCHBLK/JBLOCK
     :      /WHERE/IMCDF,NREC
*
*   This parameter determines the maximum matrix order for which the
*   LAPACK eigensolver is used
*
      PARAMETER (IOLPCK = 800)
*
*   This PARAMETER is the absolute error tolerance for the eigenvalues
*
      PARAMETER (ABSTOL = 1.0D-10)
*
*   This PARAMETER is the maximum number of doubleprecision words
*   in core due to the matrix. 8*1048576 words is 64 megabyte
*
       PARAMETER (NINCOR = 8388608)
c       PARAMETER (NINCOR = 1)
*       
*
      WRITE (24,*)
*
      IF (NCF .EQ. 1) THEN
*
*   Report
*
         WRITE (24,*) 'Trivial eigenvalue problem.'
*
*   Matrix of order 1: the trivial case; we assume that the value
*   of EAV is available
*
         CALL ALLOC (PNEVAL,1,8)
         CALL ALLOC (PNEVEC,1,8)
         EVAL(1) = 0.0D 0
         EVEC(1) = 1.0D 00
*
      ELSE
*
*   Matrix of order greater than 1; how many elements in a triangle?
*
         NDENSE = (NCF*(NCF+1))/2
*
         IF (NCF .LE. IOLPCK) THEN
*
*   Report selection of diagonaliser
*
            WRITE (24,*) 'LAPACK routine DSPEVX'
     :                 //' selected for eigenvalue problem.'
*
*   Allocate storage for the dense representation of the matrix
*
            CALL ALLOC (PNTEMT,NDENSE,8)
*
*   Initialise EMT
*
            CALL DINIT (NDENSE,0.0D 00,EMT,1)
*
*   Read the matrix into position from the disc file
*
            CALL POSNFL (IMCDF,NREC)
            CALL ALLOC (PNWORK,NCF,8)
            CALL ALLOC (PNIROW,NCF,8)
            IOFSET = 0
            DO 2 I = 1,NCF
               READ (IMCDF) NELC,ELSTO,(WORK(IR),IR = 1,NELC),
     :                                 (IROW(IR),IR = 1,NELC)
               EMT(IOFSET+IROW(1)-I+1) = WORK(1)-EAV
               DO 1 IR = 2,NELC
                  EMT(IOFSET+IROW(IR)-I+1) = WORK(IR)
    1          CONTINUE
               IOFSET = IOFSET+NCF-I+1
    2       CONTINUE
            CALL DALLOC (PNWORK)
            CALL DALLOC (PNIROW)
*
*   Find the eigenpairs
*
            NVECMN = NCF
            DO 3 I = 1,NVEC
               NVECMN = MIN (NVECMN,IVEC(I))
    3       CONTINUE
            NVEX = NVECMX-NVECMN+1
            CALL ALLOC (PONTRW,    NVEX,8)
            CALL ALLOC (PONTRZ,NCF*NVEX,8)
            CALL ALLOC (PNWORK,NCF*8   ,8)
            CALL ALLOC (PIWORK,NCF*5   ,4)
            CALL ALLOC (PIFAIL,    NVEX,4)
            CALL DSPEVX ('V','I','L',NCF,EMT,DUMMY,DUMMY,
     :                           NVECMN,NVECMX,ABSTOL,M,W,Z,NCF,
     :                           WORK,IWORK,IFAIL,INFO)
            IF (INFO .NE. 0) THEN
               PRINT *, 'MANEIG: Failure in DSPEVX [LAPACK]'
               PRINT *, ' with INFO = ',INFO,'.'
               STOP
            ENDIF
            CALL DALLOC (PNWORK)
            CALL DALLOC (PIWORK)
            CALL DALLOC (PIFAIL)
            CALL DALLOC (PNTEMT)
*
*   Store the eigenpairs in their proper positions
*
            CALL ALLOC (PNEVAL,    NVEC,8)
            CALL ALLOC (PNEVEC,NCF*NVEC,8)
            DO 4 I = 1,NVEC
               LOC = IVEC(I)
               EVAL(I) = W(LOC-NVECMN+1)
               IOFSET = NCF*(I-1)
               LOC = NCF*(LOC-NVECMN)
               CALL DCOPY (NCF,Z(LOC+1),1,EVEC(IOFSET+1),1)
    4       CONTINUE
            CALL DALLOC (PONTRW)
            CALL DALLOC (PONTRZ)
*
         ELSE
*
*   Report
*
            WRITE (24,*) 'DVDSON routine selected for'
     :                 //' for eigenvalue problem;'
*
*   Sparse or dense matrix multiply? On disc or in core?
*
            NBRKEV =  DBLE ((NCF+1)*(NCF+1))
     :               / 3.0D 00
            IF (NELMNT .LT. NBRKEV) THEN
               SPARSE = .TRUE.
               NSTORE = NELMNT+NELMNT/2+(NCF+1)/2
            ELSE
               SPARSE = .FALSE.
               NSTORE = NDENSE
            ENDIF
*
            CALL ALLOC (PNDIAG,NCF,8)
*
            IF (NSTORE .GT. NINCOR) THEN
*
*   Report
*
               WRITE (24,*) ' matrix stored on disc;'
*
*   Disk storage; necessarily sparse; one column of the matrix in
*   memory
*
               LDISC = .TRUE.
               SPARSE = .TRUE.
               CALL ALLOC (PNTEMT,NCF,8)
               CALL ALLOC (PNIROW,NCF,4)
               IMV = 1
*
*   Load diagonal 
*
               CALL POSNFL (IMCDF,NREC)
               DO 5 I = 1,NCF
                  READ (IMCDF) NELC,ELSTO,EMT(1)
                  DIAG(I) = EMT(1) - EAV
    5          CONTINUE
*
            ELSE
*
*   Core storage; load matrix into memory
*
               LDISC = .FALSE.
               CALL POSNFL (IMCDF,NREC)
               IF (SPARSE) THEN
*
*   Report
*
                  WRITE (24,*) ' matrix stored in sparse'
     :                       //' representation in core;'
*
                  IMV = 2
                  CALL ALLOC (PNTEMT,NELMNT,8)
                  CALL ALLOC (PNIROW,NELMNT,4)
                  CALL ALLOC (PIENDC,NCF+1,4)
                  IOFSET = 0
                  IENDC(0) = 0
                  DO 6 I = 1,NCF
                     READ (IMCDF) NELC,ELSTO,
     :                            (EMT(IR+IOFSET),IR = 1,NELC),
     :                           (IROW(IR+IOFSET),IR = 1,NELC)
                     DIAG(I)       = EMT(1+IOFSET) - EAV
		     EMT(1+IOFSET) = DIAG(I)
                     IOFSET = IOFSET+NELC
                     IENDC(I) = IOFSET
    6             CONTINUE
               ELSE
*
*   Report
*
                  WRITE (24,*) ' matrix stored in full'
     :                       //' representation in core;'
*
                  IMV = 3
                  CALL ALLOC (PNTEMT,NDENSE,8)
                  CALL DINIT (NDENSE,0.0D 00,EMT,1)
                  CALL ALLOC (PNWORK,NCF,8)
                  CALL ALLOC (PNIROW,NCF,4)
                  IOFSET = 0
                  DO 8 I = 1,NCF
                     READ (IMCDF) NELC,ELSTO,
     :                            (WORK(IR),IR = 1,NELC),
     :                            (IROW(IR),IR = 1,NELC)
                     DIAG(I) = WORK(1) -EAV
		     WORK(1) = DIAG(I)
                     DO 7 IR = 1,NELC
                        EMT(IOFSET+IROW(IR)-I+1) = WORK(IR)
    7                CONTINUE
                     IOFSET = IOFSET+NCF-I+1
    8             CONTINUE
                  CALL DALLOC (PNWORK)
                  CALL DALLOC (PNIROW)
               ENDIF
            ENDIF
*
*   Analyse the Hamiltonian matrix to determine its block structure
*
            CALL FNDBLK (LPRINT)
*
*   Allocate storage for workspace; see the header of DVDSON for
*   the expression below; the value of LIM can be reduced to NVECMX
*   plus a smaller number if storage is severely constrained
*
            LIM = MIN (NCF,NVECMX+20)
c           LWORK = 2*NCF*LIM+LIM*LIM+(NVECMX+10)*LIM+NVECMX
            LWORK = 2*NCF*LIM+LIM*LIM*2+11*LIM+NVECMX
            CALL ALLOC (PNWORK,LWORK,8)
            LIWORK = 6*LIM+NVECMX
            CALL ALLOC (PIWORK,LIWORK,4)
**changed by Misha 02/12/97
         CRITE = 1.0D-17
         CRITC = 1.0D-08
         CRITR = 1.0D-08
         ORTHO = max(1D-8,CRITR)
* end of changes

            MAXITR = MAX (NVECMX*100,NCF/10)
            CALL ALLOC (PJWORK,LIM,4)
*
*   Find MIN (NVECMX,NELBLK(I)) lowest eigenpairs in each block I
*
            CALL ALLOC (PNEVAL,    NVECMX,8)
            CALL ALLOC (PNEVEC,NCF*NVECMX,8)
*
            DMUNGO = 10.0D 99
            CALL DINIT (NVECMX,DMUNGO,EVAL,1)
*
            CALL ALLOC (PJWORK,LIM,4)
*
*   Compute the eigenpairs in each block
*
            DO 17 JBLOCK = 1,NBLOCK
*
               NVEX = MIN (NVECMX,NELBLK(JBLOCK))
               IF (LDISC) THEN
                  MBLOCK = NVEX
               ELSE
                  MBLOCK = 1
               ENDIF
               NEND = NCF*NVEX
*
               ILOW = 1
               IHIGH = NVEX
               NIV = NVEX
*
*   Initial estimates for eigenvectors
*
*changed by Misha 02/12/97

*
*   Call Davidson eigensolver
*
        IF     (IMV .EQ. 1) THEN
********************* sparse and matrix on disk **********************
	    call INIESTSD(NCF,NELBLK(JBLOCK),NIV,EMT,IENDC(0),IROW,
     :                WORK, IBLOCK, JBLOCK,IMCDF,NREC, EAV)

               CALL GDVD (SPODMV,NCF,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)
     
        ELSEIF (IMV .EQ. 2) THEN
********************* sparse and matrix in memory ********************
	    call INIEST(NCF,NELBLK(JBLOCK),NIV,EMT,IENDC(0),IROW,WORK,
     :                 IBLOCK, JBLOCK)

              CALL GDVD (SPICMV,NCF,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)


        ELSEIF (IMV .EQ. 3) THEN
**************************** dense and in memory **********************
	    call INIESTDM(NCF,NELBLK(JBLOCK),NIV,EMT,WORK,
     :                 IBLOCK, JBLOCK)

              CALL GDVD (DNICMV,NCF,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)


               ENDIF
*
               CALL CONVRT (NLOOPS,CNUM,LENTH)
               WRITE (24,*) ' '//CNUM(1:LENTH)//' iterations;'
               CALL CONVRT (NMV,CNUM,LENTH)
               WRITE (24,*) ' '//CNUM(1:LENTH)
     :                    //' matrix-vector multiplies.'
*
               IF (IERR .NE. 0) THEN
                  PRINT *, 'MANEIG: Returned from DVDSON with'
                  PRINT *, ' IERR = ',IERR,'.'
                  STOP
               ENDIF
*
*   Put the eigenpairs in order, overwriting as necessary
*
               IF (JBLOCK .EQ. 1) THEN
                  CALL DCOPY (    NVEX,WORK(NEND+1),1,
     :                                EVAL,1)
                  CALL DCOPY (NCF*NVEX,WORK(     1),1,
     :                                   EVEC,1)
               ELSE
                  DO 16 I = 1,NVEX
                     EVALI = WORK(NEND+I)
                     DO 15 J = I,NVECMX
                        IF (EVALI .LT. EVAL(J)) THEN
                           DO 14 K = NVECMX,J+1,-1
                              EVAL(K) = EVAL(K-1)
                              CALL DCOPY (NCF,
     :                                EVEC((K-2)*NCF+1),1,
     :                                EVEC((K-1)*NCF+1),1)
   14                      CONTINUE
                           EVAL(J) = EVALI
                           CALL DCOPY (NCF,
     :                               WORK((I-1)*NCF+1),1,
     :                               EVEC((J-1)*NCF+1),1)
                           GOTO 16
                        ENDIF
   15                CONTINUE
   16             CONTINUE
               ENDIF
*
   17       CONTINUE
*
            CALL DALLOC (PNTEMT)
            CALL DALLOC (PNDIAG)
            CALL DALLOC (PNWORK)
            CALL DALLOC (PIWORK)
            CALL DALLOC (PJWORK)
*
*   Rearrange and reallocate storage for the eigenpairs
*   as necessary
*
	    IF (NVEC .LT. NVECMX) THEN
               CALL ALLOC (PIWORK,NVECMX,4)
               DO 18 I = 1,NVECMX
                  IWORK(I) = I
   18          CONTINUE
               DO 19 I = 1,NVEC
                  IOFSET = IVEC(I)
                  LOC = IWORK(I)
                  IF (IOFSET .NE. LOC) THEN
                     CALL DSWAP (1,EVAL(IOFSET),1,
     :                             EVAL(I     ),1)
                     IWORK(I) = IWORK(IOFSET)
                     IWORK(IOFSET) = LOC
                     IOFSET = NCF*(IOFSET-1)
                     LOC = NCF*(I-1)
                     CALL DSWAP (NCF,EVEC(IOFSET+1),1,
     :                               EVEC(LOC   +1),1)
                  ENDIF
   19          CONTINUE
               CALL DALLOC (PIWORK)
               CALL RALLOC (PNEVAL,    NVECMX,    NVEC,8)
               CALL RALLOC (PNEVEC,NCF*NVECMX,NCF*NVEC,8)
	    END IF
*
            IF (SPARSE) THEN
               CALL DALLOC (PNIROW)
               IF (.NOT. LDISC) CALL DALLOC (PIENDC)
            ENDIF
            CALL DALLOC (PIBLOC)
            CALL DALLOC (PNELBL)
*
         ENDIF
*
      ENDIF
*
*   Allocate storage for eigenvector symmetry arrays
*
      CALL ALLOC (PIATJP,NVEC,4)
      CALL ALLOC (PIASPA,NVEC,4)
*
*   Clean up eigenvectors; determine their J/P values
*
      DO 23 J = 1,NVEC
*
*   Find the dominant component of each eigenvector
*
         IOFSET = (J-1)*NCF
*
         AMAX = 0.0D 00
         DO 20 I = 1,NCF
            WA = ABS (EVEC(I+IOFSET))
            IF (WA .GT. AMAX) THEN
               AMAX = WA
               IA = I
            ENDIF
   20    CONTINUE
*
*   Find the angular momentum and parity of the dominant component
*
         ITJPIA = ITJPO (IA)
         ISPAIA = ISPAR (IA)
*
*   Assign these to the ASF
*
         IATJPO(J) = ITJPIA
         IASPAR(J) = ISPAIA
*
*   Change sign of eigenvactor if dominant component is negative
*
         IF (EVEC(IA+IOFSET) .LT. 0.0D 00) THEN
*
*   Perform replacement with normalization and/or inversion
*
           DO 22 I = 1,NCF
              EVEC(I+IOFSET) =  -EVEC(I+IOFSET)
   22      CONTINUE
         ENDIF
*
   23 CONTINUE
*
      RETURN
      END
