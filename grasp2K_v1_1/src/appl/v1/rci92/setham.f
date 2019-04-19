************************************************************************
*                                                                      *
      SUBROUTINE SETHAM (ELSTO,ICSTRT)
*                                                                      *
*   Sets up the Hamiltonian matrix and determines the average energy.  *
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
*               [RCI92]: BRINT, IABINT, KEINT, RKINTC, VINT, VPINT.    *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
C     :        PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
C     :        PCOEIL,PCOEVL,
C     :        PCTEIL,PCTEVL,
C     :        PNEVAL,
C     :        PINDKE,PVALKE,
C     :        PIENDC,
C     :        PNTRIQ,
C     :        PNTJQS,PNJCUP,
C     :        PVINIL,PVINVL,
C     :        PINDVP,PVALVP,
C     :        PNTRPF,PNTRQF
      POINTER (PINDT1,INDT1DUMMY)
      POINTER (PINDT2,INDT2DUMMY)
      POINTER (PINDT3,INDT3DUMMY)
      POINTER (PINDT4,INDT4DUMMY)
      POINTER (PINDT5,INDT5DUMMY)
      POINTER (PINDT6,INDT6DUMMY)
      POINTER (PVALT1,VALT1DUMMY)
      POINTER (PVALT2,VALT2DUMMY)
      POINTER (PVALT3,VALT3DUMMY)
      POINTER (PVALT4,VALT4DUMMY)
      POINTER (PVALT5,VALT5DUMMY)
      POINTER (PVALT6,VALT6DUMMY)
      POINTER (PCOEIL,COEILDUMMY)
      POINTER (PCOEVL,COEVLDUMMY)
      POINTER (PCTEIL,CTEILDUMMY)
      POINTER (PCTEVL,CTEVLDUMMY)
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PINDKE,INDKEDUMMY)
      POINTER (PVALKE,VALKEDUMMY)
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNJCUP,JCUPDUMMY)
      POINTER (PVINIL,VINILDUMMY)
      POINTER (PVINVL,VINVLDUMMY)
      POINTER (PINDVP,INDVPDUMMY)
      POINTER (PVALVP,VALVPDUMMY)
      POINTER (PNTRPF,RPFDUMMY)
      POINTER (PNTRQF,RQFDUMMY)
      POINTER (PNIVEC,NVECMXDUMMY)

      POINTER (PCTEVLRK,VALTEIRK(1))                                  
      POINTER (PCTEILRK, INDTEIRK(1))
      POINTER (PNEVEC1,EVEC1(1))

      LOGICAL FIRST,FRSTCO,FRSTCT,FRSTKI,FRSTVI,FRSTVP,
     :        LDBPA,
     :        LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES,GETYN
      CHARACTER*8 CNUM1,CNUM2
      CHARACTER*2 NH
*
      EXTERNAL BREIT,BREID,COR,CORD
*
      DIMENSION TSHELL(NNNW)
*
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
      POINTER (PNTEMT,EMT(1))
      POINTER (PNIROW,IROW(1))
*
      COMMON/BCORE/ICORE(NNNW)
     :      /BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /COEILS/NDCOEA,NCOEI,PCOEIL,PCOEVL,FRSTCO
     :      /CTEILS/NDCTEA,NCTEI,PCTEIL,PCTEVL,FRSTCT
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEBUGA/LDBPA(5)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /FOPARM/ICCUT
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /KEILST/NDKEA,NKEI,PINDKE,PVALKE,FRSTKI
     :      /NCDIST/ZDIST(NNNP)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /STAT/PNTJQS,PNJCUP
     :      /STOR/KEEP(2,2)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /VINLST/NDVIN,NVINTI,PVINIL,PVINVL,FRSTVI
     :      /VPILST/NDVPA,NVPI,PINDVP,PVALVP,FRSTVP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /CTEILSRK/PCTEILRK,PCTEVLRK
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /EIGVEC1/PNEVEC1
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
      IF (IPRERUN.EQ.2) THEN
      DO IPI = 1,NVEC                         
        DO IPJ = 1,NCF
          WRITE(*,*) IPI,IPJ,EVEC1(IPJ+(IPI-1)*NCF)
        ENDDO
      ENDDO
      ENDIF

*
*   Allocate storage to arrays in COMMON/BUFFER/; these are
*   used for the Coulomb and transverse two-electron integrals
*
      CALL ALCBUF (1)
*
      INC1 = 1
      INC2 = 1
*
*   Initialisations for contributions from the Dirac-Coulomb
*   operator
*
      KT = 0
      IPT = 1
*
      INCOR = 1
*
      FRSTCO = .TRUE.
      NCOEC = 0
      NCOEI = 0
*
      FRSTCT = .TRUE.
      NCTEC = 0
      NCTEI = 0
*
      IF (LTRANS) THEN
*   
*   Find 2*JMAX
*   
         J2MAX = NKJ(1)
         DO I = 2,NW
           IF (NKJ(I).GT.J2MAX) J2MAX = NKJ(I)
         ENDDO
*
*   This is the maximum numbers of orbtitals allowed in brint.f
* 
         NWMAX = NNNW
         IF (J2MAX.EQ.11) NWMAX = 114
         IF (J2MAX.EQ.12) NWMAX = 112
         IF (J2MAX.EQ.13) NWMAX = 110
         IF (J2MAX.EQ.14) NWMAX = 108
         IF (J2MAX.EQ.15) NWMAX = 106
         IF (J2MAX.EQ.16) NWMAX = 105
         IF (J2MAX.EQ.17) NWMAX = 103
         IF (J2MAX.EQ.18) NWMAX = 101
         IF (J2MAX.EQ.19) NWMAX = 100
         IF (J2MAX.GT.20) NWMAX =  90
         IF (NWMAX.LT.NW) THEN
           WRITE(*,*) 'In setham. The number of orbitals is too'
           WRITE(*,*) 'large for the brint routine'
           STOP
         ENDIF
*
*   Initialisations for transverse interaction correction
*
         DO 2 I = 1,NW
            ICORE(I) = 0
            DO 1 J = 1,NCF
               IF (ICHOP (I,J) .LE. 0) GOTO 2
    1       CONTINUE
            ICORE(I) = 1
    2    CONTINUE
*
         DO 3 I = 1,6
            FIRST(I) = .TRUE.
            NTPI(I) = 0
    3    CONTINUE
*
         NMCBP = 0
*
         NCORE = 0
      ENDIF
*
*   Initialisations for the vacuum polarisation corrections
*
      IF (LVP) THEN
*
         CALL NCHARG
         CALL VACPOL
         DO 4 K = 2,N
            ZDIST(K) = TB(K)*RP(K)
    4    CONTINUE
         FRSTVP = .TRUE.
         NVPI = 0
*
      ENDIF
*
*   Initialisations for nuclear translational energy corrections
*
      IF (EMN .GT. 0.0D 00) THEN
         ATWINV = 1.0D 00/EMN
         IF (LNMS) THEN
            FRSTKI = .TRUE.
            NKEI = 0
         ENDIF
         IF (LSMS) THEN
            FRSTVI = .TRUE.
            NVINTI = 0
         ENDIF
      ELSE
         LNMS = .FALSE.
         LSMS = .FALSE.
      ENDIF
*
*   Generate the columns of the Hamiltonian matrix
*
      DO 10 IC = ICSTRT,NCF
*
*   Initialise the counter for the number of elements in this
*   column with magnitude greater than CUTOFF
*
         NELC = 0
*
*   Generate the part of the row on and below the diagonal
*
         IR = IC
*
*   Initialise the accumulator
*
    5    ELEMNT = 0.0D 00
*
         IF (LDBPA(2)) THEN
            WRITE (99,*)
            WRITE (99,*) '                       T  (ab)'
            WRITE (99,*) '  r  s   a    b         rs'
            IBUG1 = 1
         ENDIF
*
*   Generate the integral list for the matrix element of the
*   one-body operators
*
C         WRITE(*,*) '1',IR,IC
         IF (IPRERUN.EQ.1) THEN
           INC1 = 0
           INC2 = 0
           IF (IC.LE.NCSFPRE.OR.IC.EQ.IR) THEN 
             INC1 = 1
           ENDIF
         ENDIF

         IF (IPRERUN.EQ.2) THEN
*
*   Diagonal elements are always included
*   Off diagonal elements are included only if the
*   products of the weights from the prerun are larger
*   than the cutoff.
*
           IF (IC.EQ.IR) THEN
             INC1 = 1
             INC2 = 1
           ELSE
             INC1 = 0
             INC2 = 0
           ENDIF
           DO IPI = 1,NVEC
             PRECOEFF = 
     :       DABS(EVEC1(IC+(IPI-1)*NCF)*EVEC1(IR+(IPI-1)*NCF))
C             WRITE(*,*) IPI,IC,IR,COEFFCUT1,COEFFCUT2,PRECOEFF
             IF (PRECOEFF.GT.COEFFCUT1) INC1 = 1 
             IF (PRECOEFF.GT.COEFFCUT2) INC2 = 1 
           ENDDO
         ENDIF

         IF (INC1.EQ.1) THEN
C           WRITE(*,*) '2',IR,IC
           CALL TNSRJJ (KT,IPT,IC,IR,IA,IB,TSHELL)
*
*   Accumulate the contribution from the one-body operators:
*   kinetic energy, electron-nucleus interaction; update the
*   angular integral counter
*
           IF (IA .NE. 0) THEN
             IF (IA .EQ. IB) THEN
               DO 6 IA = 1,NW
                  TCOEFF = TSHELL(IA)
                  IF (LDBPA(2)) WRITE (99,300)
     :               IC,IR,NP(IA),NH(IA),NP(IA),NH(IA),TCOEFF
                  IF (ABS (TCOEFF) .GT. CUTOFF) THEN
                     NCOEC = NCOEC+1
                     CALL IABINT (IA,IA,TEGRAL)
                     ELEMNT = ELEMNT+TEGRAL*TCOEFF
                     IF (LNMS) THEN
                        CALL KEINT (IA,IA,TEGRAL)
                        ELEMNT = ELEMNT+TEGRAL*ATWINV*TCOEFF
                     ENDIF
                     IF (LVP) THEN
                        CALL VPINT (IA,IA,TEGRAL)
                        ELEMNT = ELEMNT+TEGRAL*TCOEFF
                     ENDIF
                  ENDIF
    6          CONTINUE
            ELSE
               TCOEFF = TSHELL(1)
               IF (LDBPA(2)) WRITE (99,300)
     :            IC,IR,NP(IA),NH(IA),NP(IB),NH(IB),TCOEFF
               IF (ABS (TCOEFF) .GT. CUTOFF) THEN
                  NCOEC = NCOEC+1
                  CALL IABINT (IA,IB,TEGRAL)
                  ELEMNT = ELEMNT+TEGRAL*TCOEFF
                  IF (LNMS) THEN
                     CALL KEINT (IA,IB,TEGRAL)
                     ELEMNT = ELEMNT+TEGRAL*ATWINV*TCOEFF
                  ENDIF
                  IF (LVP) THEN
                     CALL VPINT (IA,IB,TEGRAL)
                     ELEMNT = ELEMNT+TEGRAL*TCOEFF
                  ENDIF
               ENDIF
             ENDIF
           ENDIF
*
           IBUG1 = 0
*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*
           IF (LDBPA(3)) THEN
              WRITE (99,*)
              WRITE (99,*) '                                  k'
              WRITE (99,*) '                                 V  (abcd)'
              WRITE (99,*) 'r  s   a    c    b    d   k       rs'
              IBUG1 = 1
           ENDIF
*
           NVCOEF = 0
*
           CALL RKCO (IC,IR,COR,CORD,INCOR)
*
           DO 7 I = 1,NVCOEF
             VCOEFF = COEFF(I)
             IF (ABS (VCOEFF) .GT. CUTOFF) THEN
               NCTEC = NCTEC+1
               IF (LSMS) THEN
                  IF (LABEL(5,I) .EQ. 1) THEN
                     CALL VINT (LABEL(1,I),LABEL(3,I),TGRL1)
                     CALL VINT (LABEL(2,I),LABEL(4,I),TGRL2)
                     ELEMNT = ELEMNT-TGRL1*TGRL2*ATWINV*VCOEFF
                  ENDIF
               ENDIF
               CALL RKINTC (LABEL(1,I),LABEL(2,I),
     :                      LABEL(3,I),LABEL(4,I),
     :                      LABEL(5,I),TEGRAL)
               ELEMNT = ELEMNT+TEGRAL*VCOEFF
             ENDIF
    7      CONTINUE
*
           IBUG1 = 0
         ENDIF
*
         IF (LTRANS) THEN
           IF (INC2.EQ.1) THEN
C           WRITE(*,*) '3',IR,IC
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
           IF (LDBPA(4)) THEN
             WRITE (99,*)
             WRITE (99,*)
     :          '                                     kt'
             WRITE (99,*)
     :          '                                    V  (abcd)'
             WRITE (99,*)
     :          '  r  s   a    b    c    d   k  t     rs'
             IBUG1 = 1
           ENDIF
*
           NVCOEF = 0
*
            CALL RKCO (IC,IR,BREIT,BREID,1)
*
            DO 8 I = 1,NVCOEF
               IF (ABS (COEFF(I)) .GT. CUTOFF) THEN
                  NMCBP = NMCBP+1
                  ITYPE = ABS (LABEL(6,I))
                  IF (ITYPE.EQ.1) THEN
                    CALL BRINT1 (LABEL(1,I),LABEL(2,I),
     :                              LABEL(3,I),LABEL(4,I),
     :                              LABEL(5,I),TEGRAL)
                  ELSEIF (ITYPE.EQ.2) THEN
                    CALL BRINT2 (LABEL(1,I),LABEL(2,I), 
     :                              LABEL(3,I),LABEL(4,I),
     :                              LABEL(5,I),TEGRAL)
                  ELSEIF (ITYPE.EQ.3) THEN
                    CALL BRINT3 (LABEL(1,I),LABEL(2,I),
     :                              LABEL(3,I),LABEL(4,I),
     :                              LABEL(5,I),TEGRAL)
                  ELSEIF (ITYPE.EQ.4) THEN
                    CALL BRINT4 (LABEL(1,I),LABEL(2,I),
     :                              LABEL(3,I),LABEL(4,I),
     :                              LABEL(5,I),TEGRAL)
                  ELSEIF (ITYPE.EQ.5) THEN
                    CALL BRINT5 (LABEL(1,I),LABEL(2,I),
     :                              LABEL(3,I),LABEL(4,I),
     :                              LABEL(5,I),TEGRAL)
                  ELSEIF (ITYPE.EQ.6) THEN
                    CALL BRINT6 (LABEL(1,I),LABEL(2,I),
     :                              LABEL(3,I),LABEL(4,I),
     :                              LABEL(5,I),TEGRAL)
                  ENDIF 
                  CONTR = COEFF(I)*TEGRAL
                  IF (LABEL(6,I) .GT. 0) THEN
                     ELEMNT = ELEMNT+CONTR
                  ELSE
                     NCORE = NCORE+1
                     ELSTO = ELSTO+CONTR
                  ENDIF
               ENDIF
    8       CONTINUE
*
            IBUG1 = 0
*
            IF (IR .EQ. IC) ELEMNT = ELEMNT+ELSTO
*
         ENDIF
         ENDIF
*
*   Is the contribution above threshold? If so, store it
*

         IF ( (IR.EQ.IC) .OR. (ABS (ELEMNT) .GT. CUTOFF) ) THEN
            NELC = NELC+1
            EMT(NELC) = ELEMNT
            IROW(NELC) = IR
         ENDIF
*
         IF (IR .LT. NCF) THEN
            IF (LFORDR .AND. (IC .GT. ICCUT)) GOTO 9
            IR = IR+1
            GOTO 5
         ENDIF
*
*   This column is done; write it to disk
*
    9    WRITE (26) NELC,ELSTO,(EMT(IR),IR = 1,NELC),
     :                         (IROW(IR),IR = 1,NELC)
*
         EAV = EAV + EMT(1)
*
         CALL CONVRT (IC,CNUM1,LENTH1)
         CALL CONVRT (NELC,CNUM2,LENTH2)
         PRINT *,'Column '//CNUM1(1:LENTH1)//': '
     :                    //CNUM2(1:LENTH2)//' nonzero elements;'
*
*   Update the counter for the total number of elements
*
         NELMNT = NELMNT+NELC
*
   10 CONTINUE
*
*   Deallocate storage for the arrays in /BUFFER/
*
      CALL ALCBUF (3)
*
*   Deallocate storage for the integral lists from the
*   Dirac-Coulomb operator; the storage was allocated
*   in IABINT and RKINTC
*
      IF (NCOEI .GT. 0) THEN
         CALL DALLOC (PCOEIL)
         CALL DALLOC (PCOEVL)
      ENDIF
*
Cww      IF (NCTEI .GT. 0) THEN
C         CALL DALLOC (PCTEIL)
C         CALL DALLOC (PCTEVL)
C      ENDIF
      CALL DALLOC (PCTEVLRK)
      CALL DALLOC (PCTEILRK)
*
      IF (LTRANS) THEN
*
*   Deallocate storage for the integral lists from the
*   transverse photon interaction operator; this storage
*   was allocated in BRINT
*
         IF (NTPI(1) .GT. 0) THEN
            CALL DALLOC (PINDT1)
            CALL DALLOC (PVALT1)
         ENDIF
         IF (NTPI(2) .GT. 0) THEN
            CALL DALLOC (PINDT2)
            CALL DALLOC (PVALT2)
         ENDIF
         IF (NTPI(3) .GT. 0) THEN
            CALL DALLOC (PINDT3)
            CALL DALLOC (PVALT3)
         ENDIF
         IF (NTPI(4) .GT. 0) THEN
            CALL DALLOC (PINDT4)
            CALL DALLOC (PVALT4)
         ENDIF
         IF (NTPI(5) .GT. 0) THEN
            CALL DALLOC (PINDT5)
            CALL DALLOC (PVALT5)
         ENDIF
         IF (NTPI(6) .GT. 0) THEN
            CALL DALLOC (PINDT6)
            CALL DALLOC (PVALT6)
         ENDIF
*
      ENDIF
*
*   Deallocate storage for the nuclear motional energy integral
*   lists; this was allocated in KEINT and VINT
*
      IF (LNMS) THEN
         IF (NKEI .GT. 0) THEN
            CALL DALLOC (PINDKE)
            CALL DALLOC (PVALKE)
         ENDIF
      ENDIF
      IF (LSMS) THEN
         IF (NVINTI .GT. 0) THEN
            CALL DALLOC (PVINIL)
            CALL DALLOC (PVINVL)
         ENDIF
      ENDIF
*
*   Deallocate storage for the vacuum-polarisation integral list;
*   this was allocated in VPINT
*
      IF (LVP) THEN
         IF (NVPI .GT. 0) THEN
            CALL DALLOC (PINDVP)
            CALL DALLOC (PVALVP)
         ENDIF
      ENDIF
*
*   Deallocate storage for the susbshell quantum number
*   array
*
      CALL DALLOC (PNTJQS)
*
*   Deallocate storage for the orbital wavefunction arrays if the
*   self-energy is not to be estimated; this storage was allocated
*   in SUBROUTINE LODRWF
*
      IF (.NOT. LSE) THEN
         CALL DALLOC (PNTRPF)
         CALL DALLOC (PNTRQF)
      ENDIF
*
*   Compute the density of the Hamiltonian matrix
*
      DENSTY =   DBLE (NELMNT)
     :         / DBLE ((NCF*(NCF+1))/2)
*
*   Printouts
*
      WRITE (24,301) CUTOFF
*
      WRITE (24,302) NCOEI
      WRITE (24,303) NCOEC
      WRITE (24,304) NCTEI
      WRITE (24,305) NCTEC
*
      IF (LTRANS) THEN
         WRITE (24,306) NTPI
         WRITE (24,307) NMCBP
         WRITE (24,308) NCORE
      ENDIF
*
      IF (LVP) WRITE (24,309) NVPI
*
      IF (LNMS) WRITE (24,310) NKEI
      IF (LSMS) WRITE (24,311) NVINTI
*
      WRITE (24,312) NELMNT
      WRITE (24,313) DENSTY

!      ...Moved to matrix.f
!      EAV = EAV / DBLE (NCF)

      RETURN
*
  300 FORMAT (2(1X,1I2),2(1X,1I2,1A2),1X,1PD19.12)
  301 FORMAT ('CUTOFF set to ',1PD22.15)
  302 FORMAT ('Number of Dirac-Coulomb one-electron radial integrals',
     :        ' computed and stored: ',1I8)
  303 FORMAT ('Number of one-electron angular integrals that exceed',
     :        ' CUTOFF: ',1I8)
  304 FORMAT ('Number of Coulomb two-electron radial integrals',
     :        ' computed and stored: ',1I8)
  305 FORMAT ('Number of two-electron angular integrals that exceed',
     :        ' CUTOFF: ',1I8)
  306 FORMAT ('Number of transverse two-electron radial integrals',
     :        ' computed and stored: ',6I8)
  307 FORMAT ('Number of MCBP coefficients that exceed',
     :        ' CUTOFF: ',1I8)
  308 FORMAT ('Number of core coefficients that exceed',
     :        ' CUTOFF: ',1I8)
  309 FORMAT ('Number of vacuum polarisation integrals computed and',
     :        ' stored: ',1I8)
  310 FORMAT ('Number of kinetic energy integrals computed and',
     :        ' stored: ',1I8)
  311 FORMAT ('Number of Vinti integrals computed and',
     :        ' stored: ',1I8)
  312 FORMAT ('Number of elements that exceed CUTOFF in the lower',
     :        ' triangle of this Hamiltonian matrix: ',1I8)
  313 FORMAT ('Density of this Hamiltonian matrix: ',1PD22.15)
*
      END
