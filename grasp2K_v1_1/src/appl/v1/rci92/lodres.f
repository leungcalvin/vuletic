************************************************************************
*                                                                      *
      SUBROUTINE LODRES
*                                                                      *
*   Loads the data from the  .res  file. A number of checks are made   *
*   to ensure correctness and consistency.                             *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, GETYN, SETQIC.                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES
      CHARACTER*256 RECORD
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT

      COMMON/iounit/istdi,istdo,istde
*
*   Entry message
*
      WRITE(istde,*) 'Loading REStart File ...'
*
*   Read the basic parameters of the electron cloud; check these
*   against those deduced from the  .csl file
*
      READ (26) NELECR,NCFRES,NWRES
*
      IF (NELECR .NE. NELEC) THEN
         WRITE(istde,*) 'LODRES: Number of electrons does not'
         WRITE(istde,*) ' match that from Configuration List'
         WRITE(istde,*) ' File.'
         GOTO 2
      ENDIF
*
      IF (NCFRES .NE. NCF) THEN
         WRITE(istde,*) 'LODRES: Number of CSFs does not'
         WRITE(istde,*) ' match that from Configuration List'
         WRITE(istde,*) ' File.'
         GOTO 2
      ENDIF
*
      IF (NWRES .NE. NW) THEN
         WRITE(istde,*) 'LODRES: Number of subshells does not'
         WRITE(istde,*) ' match that from Configuration List'
         WRITE(istde,*) ' File.'
         GOTO 2
      ENDIF
*
*   Read the nuclear parameters
*
      READ (26) Z,EMN
      READ (26) NPARM,(PARM(I),I = 1,NPARM)
      READ (26) N,(ZZ(I),I = 1,N),NNUC
*
*   Error if the number of grid points is insufficient
*
      IF (N .GT. NNNP) THEN
         CALL CONVRT (N,RECORD,LENTH)
         WRITE(istde,*) 'LODRES: Number of grid points, '
     :           //RECORD(1:LENTH)//','
         WRITE(istde,*) ' exceeds dimensioned allowance, NNNP.'
         GOTO 2
      ENDIF
*
*   Read the physical effects specifications
*
      READ (26) C,LFORDR,ICCUT,LTRANS,WFACT,LVP,LNMS,LSMS
*
*   Read the remaining parameters controlling the radial grid and the
*   grid arrays
*
      NP10 = N+10
      READ (26) RNT,H,HP,
     :          (R(I),I = 1,NP10),
     :          (RP(I),I = 1,NP10),
     :          (RPOR(I),I = 1,NP10)
*
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Allocate storage for the radial wavefunction arrays
*
      CALL ALLOC (PNTRPF,NNNP*NW,8)
      CALL ALLOC (PNTRQF,NNNP*NW,8)
*
*   Read the orbital wavefunctions and the associated arrays
*
      DO 1 J = 1,NW
         READ (26) E(J),GAMA(J),PZ(J),MF(J)
         READ (26) (PF(I,J),I = 1,MF(J)),(QF(I,J),I = 1,MF(J))
    1 CONTINUE
*
      WRITE(istde,*) ' ... load complete;'
*
*   Determine if the self-energy contribution is to be estimated
*
      WRITE(istde,*) 'Estimate contributions from the self-energy?'
      YES = GETYN ()
      IF (YES) THEN
         LSE = .TRUE.
      ELSE
         LSE = .FALSE.
      ENDIF
*
      RETURN
*
    2 CLOSE (26)
      STOP
*
      END
