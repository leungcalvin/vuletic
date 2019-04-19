************************************************************************
*                                                                      *
      SUBROUTINE matFCO (K,ncf,IA,IB,imzero)
!-----------------------------------------------------------------------
!
! Generate a vector of FCO (the original IR replaced by (1:NCF)
! Common block is used to do the communication
! Purpose:
!  1) to vectorize do-loops involving calls to FCO from
!     rscf92/setcof.f (currently)
!  2) Of course to reduce calls to FCO
!  3) It also reduces calls to CLRX (totally in a program run, not
!     for per call to FCO or matFCO since matFCO always calls CLRX 
!     once while FCO may not call CLRX in some cases) since matFCO
!     will be called many times less than FCO
! Price to pay:
!   One double precision array of length NCF
! Effects on other subroutines:
!   Only rscf92/setcof.f is made to call this new routine since FCO
!   was heavily called there.
! Developed from rscf92/fco.f . The original comments are kept below
!
!XHH 1997-04-04
!
!-----------------------------------------------------------------------
*                                                                      *
*   This routine evaluates a coefficient                               *
*                                                                      *
*                                K                                     *
*                               f   (IA,IB)                            *
*                                IR                                    *
*                                                                      *
*   Here  K  is the multipolarity, IR  is the sequence number of the   *
*   configuration, and  IA  and  IB  are orbital  sequence  numbers.   *
*   ( I P Grant,  B J McKenzie,  P H Norrington, D F Mayers, and N C   *
*   Pyper,  Computer Phys Commun 21 (1980) 207-231, Eqs (6). )         *
*                                                                      *
*   Call(s) to: [LIB92]: CLRX, IQ.                                     *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
! . Two calls to IQ are physically inlined here to reduce the overhead.
! . In the inlined text, the checks for the integer range (ref. 
!    lib92/iq.f for details) are removed since FCO is called only by
!    SETCOF and SETHAM where the input arguments always fall in the
!    suitable range.
!XHH 1997.03.05 
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPA
      CHARACTER*2 NH
*
      COMMON/DEBUGA/LDBPA(5)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)

!XHH
CGG      PARAMETER (NNNWP = 30)
      INTEGER*4 iIQA
      POINTER (PNTRIQ,iIQA(NNNWP,1))
      COMMON/ORB2/NCF,NW,PNTRIQ

! New pointer and array for the vector
      LOGICAL imzero
      POINTER (pnfgco,fgco(1))
      COMMON/fg2cof/pnfgco

!
      imzero = .FALSE.

      IF (IA .EQ. IB) THEN

         inp2 = 8*MOD(ia-1,4)

         IF (K .EQ. 0) THEN
            DO ir = 1, ncf
               iqa = IBITS(iiqa((ia-1)/4+1,ir),inp2,8)
               FgCO(ir) = DBLE ((iqa *(iqa-1))/2)
            ENDDO

         ELSE
            IQF = NKJ(IA)+1
            KAPPA = NAK(IA)
            threej = CLRX (KAPPA,K,KAPPA)

            DO ir = 1, ncf
               iqa = IBITS(iiqa((ia-1)/4+1,ir),inp2,8)
               IF (IQA .EQ. IQF) THEN
                  FgCO(ir) = -0.5D0 * (threej*DBLE(IQA))**2
               ELSE
                  FgCO(ir) = 0.0D0
               ENDIF
            ENDDO

         ENDIF

      ELSE

         IF (K .EQ. 0) THEN
            inp2 = 8*MOD(ia-1,4)
            inp4 = 8*MOD(ib-1,4)
            DO ir = 1, ncf
               FgCO(ir) = DBLE(IBITS(iiqa((ia-1)/4+1,ir),inp2,8)
     &                *IBITS(iiqa((ib-1)/4+1,ir),inp4,8))
            ENDDO

         ELSE
            imzero = .TRUE.
            DO ir = 1, ncf
               FgCO(ir) = 0.0D0
            ENDDO
         ENDIF

      ENDIF

      RETURN
      END
