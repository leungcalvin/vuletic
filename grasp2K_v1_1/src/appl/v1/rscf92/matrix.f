************************************************************************
*                                                                      *
      SUBROUTINE MATRIX
*                                                                      *
*   Calls routines to form the Hamiltonian matrix and to diagonalise   *
*   it. The total angular momenta and parity of each  ASF  is found;   *
*   the eigenvectors are normalised so  that the sign of the largest   *
*   element is positive.                                               *
*                                                                      *
*   Call(s) to: [LIB92] ALLOC, DALLOC.                                 *
*               [RSCF92]: SETHAM, HMOUT.                               *
*               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL                       *
*                                                                      *
*                                         Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNIROW,PNTRIQ
      POINTER (PNIROW,IROWDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL FIRST,LDBPG
*
      POINTER (PNCMVL,CMVL(1))
*
      POINTER (PCDAMP,CDAMP(1))
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
*
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /DEBUG/LDBPG(5)
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /SYMA/PIATJP,PIASPA
*
      DATA FIRST /.TRUE./
*
*   Entry message; don't print for (E)AL calculations
*
      WRITE (*,300)
*
*   Save previous estimate of eigenvectors
*
      IF (.NOT. FIRST) THEN
         CALL ALLOC (PNCMVL,NCF*NCMIN,8)
         CALL DCOPY (NCF*NCMIN,EVEC,1,CMVL,1)
      ENDIF
*
*   Generate the Hamiltonian matrix
*
      CALL SETHAM
*
*   Print Hamiltonian matrix if debug option is set
*
      IF (LDBPG(4)) CALL HMOUT
*
*   Write out average energy
*
      WRITE (*,302) EAV
*
*   Subtract the average energy from the diagonal elements to reduce
*   the condition number of the matrix
*
      DO 1 I = 1,NCF
         IDIAG = IENDC(I-1)+1
         EMT(IDIAG) = EMT(IDIAG)-EAV
    1 CONTINUE
*
*   Compute and store eigenpairs
*
      CALL MANEIG (LDBPG(3))
*
*   Damp and Schmidt orthogonalise eigenvectors for OL calculations
*
      IF (.NOT. FIRST) THEN
*
         DO 6 J = 1,NCMIN
*
            IOFSET = (J-1)*NCF
            JOTHER = J
*
            CDAMPJ = CDAMP(J)
            IF (CDAMPJ .NE. 0.0D 00) THEN
*
               OMCDAJ = 1.0D 00-CDAMPJ
*
*   Determine the principal component
*
    2          AMAX = 0.0D 00
               DO 3 I = 1,NCF
                  EVECIJ = OMCDAJ*EVEC(I+IOFSET)
     :                    +CDAMPJ*CMVL(I+IOFSET)
                  EVEC(I+IOFSET) = EVECIJ
                  WA = ABS (EVECIJ)
                  IF (WA .GT. AMAX) THEN
                     AMAX = WA
                     IA = I
                  ENDIF
    3          CONTINUE
*
*   Find the angular momentum and parity of the dominant component
*
               ITJPIA = ITJPO (IA)
               ISPAIA = ISPAR (IA)
*
*   Remove any spurious component from the eigenvectors; compute
*   the normalization factor at the same time
*
               SUM = 0.0D 00
               DO 4 I = 1,NCF
                  IF ((ITJPO (I) .NE. ITJPIA) .OR.
     :                (ISPAR (I) .NE. ISPAIA)) THEN
                     EVEC(I+IOFSET) = 0.0D 00
                  ELSE
                     SUM = SUM+EVEC(I+IOFSET)**2
                  ENDIF
    4          CONTINUE
               DNFAC = 1.0D 00/SQRT (SUM)
*
*   Renormalize and/or invert as necessary
*
               IF (EVEC(IA+IOFSET) .LT. 0.0D 00)
     :            DNFAC = -DNFAC
               CALL DSCAL (NCF,DNFAC,EVEC(IOFSET+1),1)
*
*   Schmidt orthogonalise
*
    5          JOTHER = JOTHER-1
               IF (JOTHER .GE. 1) THEN
                  JOFSET = (JOTHER-1)*NCF
                  OVRLAP = DDOT (NCF,EVEC(IOFSET+1),1,
     :                                       EVEC(JOFSET+1),1)
                  IF (OVRLAP .NE. 0.0D 00) THEN
                     OMCDAJ = 1.0D 00
                     CDAMPJ = -OVRLAP
                     CALL DCOPY (NCF,EVEC(JOFSET+1),1,
     :                                       CMVL(IOFSET+1),1)
                     GOTO 2
                  ELSE
                     GOTO 5
                  ENDIF
               ENDIF
*
            ENDIF
*
    6    CONTINUE
*
      ENDIF
*
*   Write out the eigenpair information
*
*   Position the file
*
      REWIND (25)
      DO 7 I = 1,4
         READ (25)
    7 CONTINUE
*
*   Write ASF symmetries, eigenvalues, and eigenvectors to GRASP92
*   MIXing coefficients File
*
      WRITE (25) (IATJPO(I),IASPAR(I),I = 1,NCMIN)
      WRITE (25) EAV,(EVAL(I),I = 1,NCMIN)
      WRITE (25) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NCMIN)
*
*   Deallocate the temporary storage
*
      IF (.NOT. FIRST) THEN
         CALL DALLOC (PNCMVL)
      ELSE
         FIRST = .FALSE.
      ENDIF
*
      RETURN
*
  300 FORMAT (/'MATRIX ...')
  302 FORMAT (/' Average energy = ',1PD18.10,' Hartrees')
*
      END
