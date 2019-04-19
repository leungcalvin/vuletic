************************************************************************
*                                                                      *
      SUBROUTINE SETHAM
*                                                                      *
*   This  SUBROUTINE  sets up the  Hamiltonian matrix and determines   *
*   the average energy.                                                *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CLRX, DALLOC, RINTI, SLATER.           *
*               [RSCF92]: FCO, GCO, IQ.                                *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL SET
*
      POINTER (PCOEFF,COEFF(1))
      POINTER (PICLMN,ICLMN(1))
      POINTER (PINDEX,INDEX(1))
*
      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
*
      COMMON/EIGVAL/EAV,PNEVAL
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /MCPA/KMAXF
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
*   Allocate storage that is local to this subroutine
*
      NDIM = 1
      CALL ALLOC (PCOEFF,NDIM,8)
      CALL ALLOC (PICLMN,NDIM,4)
      CALL ALLOC (PINDEX,NDIM,4)
*
*   Other initializations
*
      CALL ALLOC (PNTEMT,NELMNT,8)
*
      NWM1 = NW-1
*
      DO 1 I = 1,NELMNT
         EMT(I) = 0.0D 00
    1 CONTINUE
*
*   Accumulate diagonal terms that do not require MCP coefficients
*
*   Piece involving I(a,a) integrals
*

      DO 3 IA = 1,NW
         SET = .FALSE.
         DO 2 IR = 1,NCF
            QA = DBLE (IQ (IA,IR))
            IF (QA .GT. 0.0D 00) THEN
               IF (.NOT. SET) THEN
                  DIAA = RINTI (IA,IA,0)
                  SET = .TRUE.
               ENDIF
               IDIAG = IENDC(IR-1)+1
               EMT(IDIAG) = EMT(IDIAG)+QA*DIAA
            ENDIF
    2    CONTINUE
    3 CONTINUE
*
*                    0
*   Piece involving F (a,a) integrals
*
      kzero = 0
      DO 5 IA = 1,NW
         SET = .FALSE.
         DO 4 IR = 1,NCF
            COEF = FCO (kzero,IR,IA,IA)
            IF (ABS (COEF) .GT. 0.0D 00) THEN
               IF (.NOT. SET) THEN
                  F0AA = SLATER (IA,IA,IA,IA,0)
                  SET = .TRUE.
               ENDIF
               IDIAG = IENDC(IR-1)+1
               EMT(IDIAG) = EMT(IDIAG)+COEF*F0AA
            ENDIF
    4    CONTINUE
    5 CONTINUE
*
*                    k
*   Piece involving F (a,a) integrals
*
      KM = 0
      K = 0
    6 K = K+2
      DO 8 IA = 1,NW
         K0 = NKJ(IA)-1
         IF (K0 .GT. KM) KM = K0
         IF (K .LE. K0) THEN
            SET = .FALSE.
            DO 7 IR = 1,NCF
               COEF = FCO (K,IR,IA,IA)
               IF (ABS (COEF) .GT. 0.0D 00) THEN
                  IF (.NOT. SET) THEN
                     FKAA = SLATER (IA,IA,IA,IA,K)
                     SET = .TRUE.
                  ENDIF
                  IDIAG = IENDC(IR-1)+1
                  EMT(IDIAG) = EMT(IDIAG)+COEF*FKAA
               ENDIF
    7       CONTINUE
         ENDIF
    8 CONTINUE
      IF (K .LT. KM) GOTO 6
*
*                    0
*   Piece involving F (a,b) integrals
*
      kzero = 0
      DO 11 IA = 1,NWM1
         IAP1 = IA+1
         DO 10 IB = IAP1,NW
            SET = .FALSE.
            DO 9 IR = 1,NCF
               COEF = FCO (kzero,IR,IA,IB)
               IF (ABS (COEF) .GT. 0.0D 00) THEN
                  IF (.NOT. SET) THEN
                     F0AB = SLATER (IA,IB,IA,IB,0)
                     SET = .TRUE.
                  ENDIF
                  IDIAG = IENDC(IR-1)+1
                  EMT(IDIAG) = EMT(IDIAG)+COEF*F0AB
               ENDIF
    9       CONTINUE
   10    CONTINUE
   11 CONTINUE
*
*                    k
*   Piece involving G (a,b) integrals
*
      KM = 0
      K = -1
   12 K = K+1
      DO 15 IA = 1,NWM1
         NKJIA = NKJ(IA)
         IAP1 = IA+1
         DO 14 IB = IAP1,NW
            NKJIB = NKJ(IB)
            SET = .FALSE.
            IF (NAK(IA)*NAK(IB) .GT. 0) THEN
               KMIN = ABS ((NKJIA-NKJIB)/2)
            ELSE
               KMIN = ABS ((NKJIA-NKJIB)/2)+1
            ENDIF
            IF (MOD (K-KMIN,2) .EQ. 0) THEN
               KMAX = (NKJIA+NKJIB)/2
               IF (KMAX .GT. KM) KM = KMAX
               IF ((K .GE. KMIN) .AND. (K .LE. KMAX)) THEN
                  DO 13 IR = 1,NCF
                     COEF = GCO (K,IR,IA,IB)
                     IF (ABS (COEF) .GT. 0.0D 00) THEN
                        IF (.NOT. SET) THEN
                           GKAB = SLATER (IA,IB,IB,IA,K)
                           SET = .TRUE.
                        ENDIF
                        IDIAG = IENDC(IR-1)+1
                        EMT(IDIAG) = EMT(IDIAG)+COEF*GKAB
                     ENDIF
   13             CONTINUE
               ENDIF
            ENDIF
   14    CONTINUE
   15 CONTINUE
      IF (K .LT. KM) GOTO 12
*
*   Accumulate one-electron terms that require MCP coefficients
*
      REWIND (31)
      READ (31)
      READ (31)
      READ (31)
*
*   Attempt to read another block of data
*
   16 READ (31,IOSTAT = IOS) LAB,NCONTR
*
      IF (IOS .EQ. 0) THEN
*
*   Read successful; decode the labels of I(ab)
*
         IA = MOD (LAB,KEY)
         IB = LAB/KEY
*
*   Compute I(ab)
*
         TEGRAL = RINTI (IA,IB,0)
*
*   Ensure that storage is adequate to read in the rest of
*   this block
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
*   Store all the contributions from this integral
*
         DO 17 I = 1,NCONTR
            LOC = INDEX(I)
            EMT(LOC) = EMT(LOC)+TEGRAL*COEFF(I)
   17    CONTINUE
*
*   Return to the start of the loop
*
         GOTO 16
*
      ENDIF
*
*   Accumulate two-electron terms that require MCP coefficients
*
      DO 20 NFILE = 32,32+KMAXF
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
         nnnblocks = 0
   18    READ (NFILE,IOSTAT = IOS) LAB,NCONTR
         nnnblocks = nnnblocks + 1
*
         IF (IOS .EQ. 0) THEN
*                                          k
*   Read successful; decode the labels of R (abcd)
*
            ID = MOD (LAB,KEY)
            LAB = LAB/KEY
            IB = MOD (LAB,KEY)
            LAB = LAB/KEY
            IC = MOD (LAB,KEY)
            IA = LAB/KEY
*
*   Compute I(ab)
*
            TEGRAL = SLATER (IA,IB,IC,ID,K)
*
*   Ensure that storage is adequate to read in the rest of
*   this block
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
*   Store all the contributions from this integral
*
            DO 19 I = 1,NCONTR
               LOC = INDEX(I)
               EMT(LOC) = EMT(LOC)+TEGRAL*COEFF(I)
   19       CONTINUE
*
*   Return to the start of the loop
*
            GOTO 18
*
         ENDIF
*
   20 CONTINUE
*
*   Deallocate storage that is local to this routine
*
      CALL DALLOC (PCOEFF)
      CALL DALLOC (PICLMN)
      CALL DALLOC (PINDEX)
*
*   Determine average energy
*
      EAV = 0.0D 00
      DO 21 IR = 1,NCF
         EAV = EAV+EMT(IENDC(IR-1)+1)
   21 CONTINUE
      EAV = EAV/DBLE (NCF)
*
      RETURN
      END
