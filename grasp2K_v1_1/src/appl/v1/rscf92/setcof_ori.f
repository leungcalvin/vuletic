************************************************************************
*                                                                      *
      SUBROUTINE SETCOF (EOL,J)
*                                                                      *
*   This  subroutine  sets  up the coefficients and orbital pointers   *
*   for the direct and exchange potentials for orbital  J .  It also   *
*   sets  up  the  coefficients  and  pointers for the inhomogeneous   *
*   terms arising from off-diagonal I (a,b) integrals.                 *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
*               [RSCF92]: ALCSCA, DSUBRS.                              *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 21 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTEMT,PIENDC,PNTRIQ
      POINTER (PNTEMT,EMTDUMMY)
      POINTER (PIENDC,IENDCDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL EOL
*
      POINTER (PCOEFF,COEFF(1))
      POINTER (PICLMN,ICLMN(1))
      POINTER (PINDEX,INDEX(1))
*
      POINTER (PNIROW,IROW(1))
      POINTER (PNTRDA,DA(1))
      POINTER (PNTRXA,XA(1))
      POINTER (PNTRYA,YA(1))
      POINTER (PNTNDA,NDA(1))
      POINTER (PNTNXA,NXA(1))
      POINTER (PNTNYA,NYA(1))
*
      DIMENSION INDEXS(4)
*
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /MCPA/KMAXF
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /SCF1/UCF(NNNW)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
*
      PARAMETER (EPS = 1.0D-10)
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
*   Initializations
*
      NDIM = 1
      CALL ALLOC (PCOEFF,NDIM,8)
      CALL ALLOC (PICLMN,NDIM,4)
      CALL ALLOC (PINDEX,NDIM,4)
*
      NDCOF = 0
      NXCOF = 0
      NYCOF = 0
*
      NAKJ = NAK(J)
      NKJJ = NKJ(J)
      UCFJ = UCF(J)
*
*   Generate YA coefficients that do not require MCP output list
*
      ILABEL = 0
      NWTERM = KEY*(KEY+1)
      DO 3 IB = 1,NW
         ILABEL = ILABEL+NWTERM
         IF (IB .EQ. J) THEN
            KMAX = NKJJ-1
         ELSE
            KMAX = 0
         ENDIF
         DO 2 K = 0,KMAX,2
            SUMR = 0.0D 00
            DO 1 IR = 1,NCF
               SUMR = SUMR+DSUBRS (EOL,IR,IR)*FCO (K,IR,J,IB)
    1       CONTINUE
            IF (IB .EQ. J) THEN
               YKAB = 2.0D 00*SUMR/UCFJ
            ELSE
               YKAB = SUMR/UCFJ
            ENDIF
            IF (ABS (YKAB) .GT. EPS) THEN
               NYCOF = NYCOF+1
               IF (NYCOF .GT. NYDIM) THEN
                  IF (NYDIM .GT. 0) THEN
                     CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,2)
                  ELSE
                     CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,1)
                  ENDIF
               ENDIF
               YA(NYCOF) = YKAB
               NYA(NYCOF) = K+ILABEL
            ENDIF
    2    CONTINUE
    3 CONTINUE
*
*   Generate XA coefficients that do not require MCP output list
*
      ILABEL = KEY*J
      NWTERM = KEY*KEY*(KEY+1)
      DO 6 IB = 1,NW
         ILABEL = ILABEL+NWTERM
         IF (IB .NE. J) THEN
            NKJIB = NKJ(IB)
            IF (NAKJ*NAK(IB) .GT. 0) THEN
               KMIN = ABS ((NKJJ-NKJIB)/2)
            ELSE
               KMIN = ABS ((NKJJ-NKJIB)/2)+1
            ENDIF
            KMAX = (NKJJ+NKJIB)/2
            DO 5 K = KMIN,KMAX,2
               SUMR = 0.0D 00
               DO 4 IR = 1,NCF
                  SUMR = SUMR+DSUBRS (EOL,IR,IR)*GCO (K,IR,J,IB)
    4          CONTINUE
               XKAB = SUMR/UCFJ
               IF (ABS (XKAB) .GT. EPS) THEN
                  NXCOF = NXCOF+1
                  IF (NXCOF .GT. NXDIM) THEN
                     IF (NXDIM .GT. 0) THEN
                        CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,2)
                     ELSE
                        CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,1)
                     ENDIF
                  ENDIF
                  XA(NXCOF) = XKAB
                  NXA(NXCOF) = K+ILABEL
               ENDIF
    5       CONTINUE
         ENDIF
    6 CONTINUE
*
*   Generate DA coefficients; these arise from the one-electron
*   integrals
*
      REWIND (31)
      READ (31)
      READ (31)
      READ (31)
*
*   Attempt to read another block of data
*
    7 READ (31,IOSTAT = IOS) LAB,NCONTR
*
      IF (IOS .EQ. 0) THEN
*
*   Read successful; decode the labels of I(ab)
*
         IA = MOD (LAB,KEY)
         IB = LAB/KEY
*
         IF ((IA .EQ. J) .OR. (IB .EQ. J)) THEN
*
*   The required subshell index has been found; ensure that storage is
*   adequate to read in the rest of this block
*
            IF (NCONTR .GT. NDIM) THEN
               CALL DALLOC (PCOEFF)
               CALL DALLOC (PICLMN)
               CALL DALLOC (PINDEX)
               NDIM = NCONTR
               CALL ALLOC (PCOEFF,NDIM,8)
               CALL ALLOC (PICLMN,NDIM,4)
               CALL ALLOC (PINDEX,NDIM,4)
            ENDIF
*
*   Read the column index, the sparse matrix index, and the
*   coefficient for all contributions from this integral
*
            READ (31) (ICLMN(I),INDEX(I),COEFF(I),I = 1,NCONTR)
*
*   Add up all the contributions from this integral; off-diagonal
*   contributions have double the weight
*
            SUM = 0.0D 00
            DO 8 I = 1,NCONTR
               IR = IROW(INDEX(I))
               IC = ICLMN(I)
               CONTR = DSUBRS (EOL,IR,IC)*COEFF(I)
               IF (IR .NE. IC) CONTR = CONTR+CONTR
               SUM = SUM+CONTR
    8       CONTINUE
*
*   Divide by twice generalized occupation number
*
            SUM = 0.5D 00*SUM/UCFJ
*
*   Put the coefficient in the list; ensure that storage is adequate
*
            NDCOF = NDCOF+1
            IF (NDCOF .GT. NDDIM) THEN
               IF (NDDIM .GT. 0) THEN
                  CALL ALCSCA (PNTNDA,PNTRDA,NDDIM,2)
               ELSE
                  CALL ALCSCA (PNTNDA,PNTRDA,NDDIM,1)
               ENDIF
            ENDIF
*
            DA(NDCOF) = SUM
            IF (IA .EQ. J) THEN
               NDA(NDCOF) = IB
            ELSE
               NDA(NDCOF) = IA
            ENDIF
*
         ELSE
*
*   No contributions from this integral; skip the rest of
*   this block of data
*
            READ (31)
*
         ENDIF
*
*   Return to the start of the loop
*
         GOTO 7
*
      ENDIF
*
*   Generate YA and XA coefficients; these arise from the two-electron
*   integrals
*
      DO 16 NFILE = 32,32+KMAXF
*
         REWIND (NFILE)
         READ (NFILE)
         READ (NFILE)
         READ (NFILE)
*
*   The multipolarity of the integral can be deduced from the file
*   unit number
*
         K = NFILE-32
*
*   Attempt to read another block of data
*
    9    READ (NFILE,IOSTAT = IOS) LAB,NCONTR
*
         IF (IOS .EQ. 0) THEN
*                                          k
*   Read successful; decode the labels of R (abcd)
*
            INDEXS(4) = MOD (LAB,KEY)
            LAB = LAB/KEY
            INDEXS(2) = MOD (LAB,KEY)
            LAB = LAB/KEY
            INDEXS(3) = MOD (LAB,KEY)
            INDEXS(1) = LAB/KEY
*
*   Determine the number of indices that match
*
            IRANK = 0
            DO 10 I = 1,4
               IF (INDEXS(I) .EQ. J) IRANK = IRANK+1
   10       CONTINUE
*
            IF (IRANK .GT. 0) THEN
*
*   At least one subshell index matches; ensure that storage
*   is adequate to read in the rest of this block
*
               IF (NCONTR .GT. NDIM) THEN
                  CALL DALLOC (PCOEFF)
                  CALL DALLOC (PICLMN)
                  CALL DALLOC (PINDEX)
                  NDIM = NCONTR
                  CALL ALLOC (PCOEFF,NDIM,8)
                  CALL ALLOC (PICLMN,NDIM,4)
                  CALL ALLOC (PINDEX,NDIM,4)
               ENDIF
*
*   Read the column index, the sparse matrix index, and the
*   coefficient for all contributions from this integral
*
               READ (NFILE) (ICLMN(I),INDEX(I),COEFF(I),I = 1,NCONTR)
*
*   Add up all the contributions from this integral; off-diagonal
*   contributions have double the weight
*
               SUM = 0.0D 00
               DO 11 I = 1,NCONTR
                  IR = IROW(INDEX(I))
                  IC = ICLMN(I)
                  CONTR = DSUBRS (EOL,IR,IC)*COEFF(I)
                  IF (IR .NE. IC) CONTR = CONTR+CONTR
                  SUM = SUM+CONTR
   11          CONTINUE
*
*   Divide by twice generalized occupation number
*
               SUM = 0.5D 00*SUM/UCFJ
*
               IF (IRANK .EQ. 1) THEN
*
*   One matching index: exchange potential contribution
*
                  NXCOF = NXCOF+1
                  IF (NXCOF .GT. NXDIM) THEN
                     IF (NXDIM .GT. 0) THEN
                        CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,2)
                     ELSE
                        CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,1)
                     ENDIF
                  ENDIF
                  XA(NXCOF) = SUM
                  DO 12 IIND = 1,4
                     IF (INDEXS(IIND) .EQ. J) THEN
                        IL = IIND+2
                        IF (IL .GT. 4) IL = IL-4
                        IORB = INDEXS(IL)
                        IL = IIND+1
                        IF (IL .GT. 4) IL = IL-4
                        IYO1 = INDEXS(IL)
                        IL = IIND+3
                        IF (IL .GT. 4) IL = IL-4
                        IYO2 = INDEXS(IL)
                        NXA(NXCOF) = (((IORB*KEY+IYO2)*KEY)+IYO1)*KEY+K
                        GOTO 9
                     ENDIF
   12             CONTINUE
               ELSEIF (IRANK .EQ. 2) THEN
*
*   Two matching indices: either direct or exchange potential
*   contribution
*
                  IFOUND = 0
                  DO 13 IIND = 1,4
                     IF (INDEXS(IIND) .EQ. J) THEN
                        IF (IFOUND .EQ. 0) THEN
                           LOC1 = IIND
                           IFOUND = IFOUND+1
                        ELSEIF (IFOUND .EQ. 1) THEN
                           LOC2 = IIND
                           GOTO 14
                        ENDIF
                     ENDIF
   13             CONTINUE
   14             IF (LOC2-LOC1 .EQ. 2) THEN
*
*   Direct contribution
*
                     NYCOF = NYCOF+1
                     IF (NYCOF .GT. NYDIM) THEN
                        IF (NYDIM .GT. 0) THEN
                           CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,2)
                        ELSE
                           CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,1)
                        ENDIF
                     ENDIF
                     YA(NYCOF) = 2.0D 00*SUM
                     IL = LOC1+3
                     IF (IL .GT. 4) IL = IL-4
                     IYO2 = INDEXS(IL)
                     IL = LOC1+1
                     IF (IL .GT. 4) IL = IL-4
                     IYO1 = INDEXS(IL)
                     NYA(NYCOF) =  (IYO2*KEY+IYO1)*KEY+K
*
                  ELSE
*
*   Exchange contribution
*
                     NXCOF = NXCOF+1
                     IF (NXCOF .GT. NXDIM) THEN
                        IF (NXDIM .GT. 0) THEN
                           CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,2)
                        ELSE
                           CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,1)
                        ENDIF
                     ENDIF
                     XA(NXCOF) = 2.0D 00*SUM
                     IL = LOC1+2
                     IF (IL .GT. 4) IL = IL-4
                     IORB = INDEXS(IL)
                     IL = LOC1+1
                     IF (IL .GT. 4) IL = IL-4
                     IYO1 = INDEXS(IL)
                     IL = LOC1+3
                     IF (IL .GT. 4) IL = IL-4
                     IYO2 = INDEXS(IL)
                     NXA(NXCOF) = (((IORB*KEY+IYO2)*KEY)+IYO1)*KEY+K
                  ENDIF
*
               ELSEIF (IRANK .EQ. 3) THEN
*
*   Three matching indices: direct and exchange potential contributions
*
                  NYCOF = NYCOF+1
                  IF (NYCOF .GT. NYDIM) THEN
                     IF (NYDIM .GT. 0) THEN
                        CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,2)
                     ELSE
                        CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,1)
                     ENDIF
                  ENDIF
                  YA(NYCOF) = 2.0D 00*SUM
                  NXCOF = NXCOF+1
                  IF (NXCOF .GT. NXDIM) THEN
                     IF (NXDIM .GT. 0) THEN
                        CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,2)
                     ELSE
                        CALL ALCSCA (PNTNXA,PNTRXA,NXDIM,1)
                     ENDIF
                  ENDIF
                  XA(NXCOF) = SUM
                  DO 15 IIND = 1,4
                     IF (INDEXS(IIND) .NE. J) THEN
                        INDIND = INDEXS(IIND)
                        IYO2 = INDIND
                        IYO1 = J
                        NYA(NYCOF) = (IYO2*KEY+IYO1)*KEY+K
                        IORB = INDIND
                        IYO1 = J
                        IYO2 = J
                        NXA(NXCOF) = (((IORB*KEY+IYO2)*KEY)+IYO1)*KEY+K
                        GOTO 9
                     ENDIF
   15             CONTINUE
*
               ELSEIF (IRANK .EQ. 4) THEN
*
*   Four matching indices: direct potential contribution
*
                  NYCOF = NYCOF+1
                  IF (NYCOF .GT. NYDIM) THEN
                     IF (NYDIM .GT. 0) THEN
                        CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,2)
                     ELSE
                        CALL ALCSCA (PNTNYA,PNTRYA,NYDIM,1)
                     ENDIF
                  ENDIF
                  YA(NYCOF) = 4.0D 00*SUM
                  IYO2 = J
                  IYO1 = J
                  NYA(NYCOF) =  (IYO2*KEY+IYO1)*KEY+K
               ENDIF
*
            ELSE
*
               READ (NFILE)
*
            ENDIF
*
            GOTO 9
*
         ENDIF
*
   16 CONTINUE
*
*   Deallocate storage for arrays local to this routine
*
      CALL DALLOC (PCOEFF)
      CALL DALLOC (PICLMN)
      CALL DALLOC (PINDEX)
*
      RETURN
      END
