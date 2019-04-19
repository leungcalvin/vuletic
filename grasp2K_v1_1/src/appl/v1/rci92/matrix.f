************************************************************************
*                                                                      *
      SUBROUTINE MATRIX
*                                                                      *
*   This SUBROUTINE calls routines to  form the  Hamiltonian  matrix   *
*   and to diagonalise it.   The total angular momenta and parity of   *
*   the ASF are found;  the eigenvectors are chosen so that the sign   *
*   of the largest element is positive.                                *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, HMOUT, IQ, POSNFL,     *
*                        WGHTD5.                                       *
*               [RCI92]: MANEIG, QED, SETHAM.                          *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PIENDC,PNTRIQ,PNTJQS,PNJCUP,PNTRPF,PNTRQF
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNJCUP,JCUPDUMMY)
      POINTER (PNTRPF,RPFDUMMY)
      POINTER (PNTRQF,RQFDUMMY)
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,LDBPG
      CHARACTER*8 CNUM
*
      DIMENSION SLFINT(NNNW)
*
      POINTER (PNETOT,ETOT(1))
*
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNEVEC1,EVEC1(1))
      POINTER (PNTEMT,EMT(1))
      POINTER (PNIROW,IROW(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
*
      COMMON/DEBUGG/LDBPG(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /STAT/PNTJQS,PNJCUP
     :      /SYMA/PIATJP,PIASPA
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WHERE/IMCDF,NREC
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /EIGVEC1/PNEVEC1
*
      WRITE (24,*)
*
*   Allocate storage for one column of the sparse representation
*   of the Hamiltonian matrix
*
      CALL ALLOC (PNTEMT,NCF,8)
      CALL ALLOC (PNIROW,NCF,4)
*
*   EAV stores the average energy
*
      EAV = 0.0D 00
*
*   Record NREC+1 stores the first column of the Hamiltonian
*   matrix
*
      NREC = 6+2*NW
*
*   Determine the number of elements that exceed CUTOFF and the
*   number of columns of the Hamiltonian matrix that have already
*   been written to the RCI92 REStart File
*
      NELMNT = 0
      ICSTRT = 0
      elsto  = 0.d0
*
      CALL POSNFL (26,NREC)
*
    1 READ (26,IOSTAT = IOS) NELC,STOEL,(EMT(IR),IR = 1,NELC),
     :                                 (IROW(IR),IR = 1,NELC)
      IF (IOS .EQ. 0) THEN
         ELSTO = STOEL
         EAV = EAV+EMT(1)
         NELMNT = NELMNT+NELC
         ICSTRT = ICSTRT+1
         GOTO 1
      ENDIF
*
      IF (ICSTRT .NE. NCF) THEN
*
*   Part of the Hamiltonian matrix remains to be generated
*
         IF (ICSTRT .GT. 0) THEN
            CALL CONVRT (ICSTRT,CNUM,LENTH)
            PRINT     *, CNUM(1:LENTH)//' columns of the'
     :                  //' Hamiltonian matrix read from'
     :                  //' the RCI92 REStart File;'
            WRITE (24,*) CNUM(1:LENTH)//' columns of the'
     :                  //' Hamiltonian matrix read from'
     :                  //' the RCI92 REStart File;'
         ENDIF
*
*   Position the RCI92 REStart File to the end of the last complete
*   record
*
         CALL POSNFL (26,NREC+ICSTRT)
*
         ICSTRT = ICSTRT+1
*
*   Generate the remainder of the Hamiltonian matrix
*
         CALL SETHAM (ELSTO,ICSTRT)
*
      ELSE
*
*   The entire Hamiltonian matrix has already been generated
*
         PRINT     *, 'Entire Hamiltonian matrix read from disc;'
         WRITE (24,*) 'Entire Hamiltonian matrix read from disc;'
*
      ENDIF
*
*   Print the Hamiltonian matrix if the debug option is set
*
      IF (LDBPG(4)) THEN
         CALL POSNFL (26,NREC)
         WRITE (99,*)
         WRITE (99,*) 'Elements of the Hamiltonian matrix:'
         WRITE (99,*)
         CALL HMOUT (26)
      ENDIF
*
*   Deallocate storage
*
      CALL DALLOC (PNTEMT)
      CALL DALLOC (PNIROW)
*
*   Allocate and deallocate memory for the mixing coefficients from the prerun
*
      IF (IPRERUN.EQ.1) CALL ALLOC (PNEVEC1,NCF*NVEC,8)
      IF (IPRERUN.EQ.2) CALL DALLOC (PNEVEC1)
*
*   Write out the average energy
*
      EAV = EAV / DBLE (NCF)
      WRITE (24,*)
      WRITE (24,300) EAV
*
*   Determine the eigenpairs
*
      IMCDF = 26
      CALL MANEIG (LDBPG(3))
*
*   Write out the eigenvalues and the dominant components of
*   the eigenvectors
*
      CALL ENGOUT (EAV,EVAL,IATJPO,IASPAR,IVEC,NVEC,3)
      CALL WGHTD5
*
*   Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
*   MIXing coefficients File; close the file; print a report
*
      WRITE (25) (IATJPO(I),IASPAR(I),I = 1,NVEC)
      WRITE (25) EAV,(EVAL(I),I = 1,NVEC)
      WRITE (25) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
*
*   Save the mixing coefficients from the prerun
*
      IF (IPRERUN.EQ.1) THEN
        DO J = 1,NVEC
          DO I = 1,NCF
            EVEC1(I+(J-1)*NCF) = EVEC(I+(J-1)*NCF) 
          ENDDO
        ENDDO
      ENDIF
*
      CLOSE (25)
*
      PRINT *, 'RCI92 MIXing coefficients File generated.'
*
*   Estimate diagonal self-energy corrections; the stored
*   eigenvalues and eigenvectors are not modified by these
*   estimates
*
      IF (LSE) THEN
         CALL ALLOC (PNETOT,NVEC,8)
         PRINT *, 'Entering QED ...'
!         CALL QED (SLFINT)
!         PRINT *, ' ... complete.'
         DO 3 J = 1,NVEC
!
! Moved here inside the eigenvalue loop
!
            CALL QED (j,SLFINT)

            ELEMNT = 0.0D 00
            IC = IVEC(J)
            DO 2 I = 1,NW
               ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
    2       CONTINUE
            ETOT(J) = EVAL(J)+ELEMNT
    3    CONTINUE
         WRITE (24,*)
         WRITE (24,*) 'Self-energy corrections estimated'
     :              //' --- these do not influence the data'
     :              //' in the RCI92 MIXing coefficients File.'
         CALL ENGOUT (EAV,ETOT,IATJPO,IASPAR,IVEC,NVEC,MODE)
         CALL DALLOC (PNETOT)
         CALL DALLOC (PNTRPF)
         CALL DALLOC (PNTRQF)
      ENDIF
*
      CALL DALLOC (PNJCUP)
      CALL DALLOC (PNTRIQ)
*
      RETURN
*
  300 FORMAT ('Average energy = ',1PD19.12,' Hartrees.')
*
      END
