************************************************************************
*                                                                      *
      SUBROUTINE gcovec (K,IA,IB)
!----------------------------------------------------------------------- 
!
!   This routine computes 
!
!      DO IR = 1, NCF
!                     K
!         fgco(IR) = g  (IA,IB)
!                     IR
!      ENDDO
!
!   Directly developed from GCO where the original comments follow
!
!XHH 1997.04.08
!
!                         !!!  _Important_ note  !!!
!
!
!     The same common block (and array) as in fcovec.f is used here
!     to store coefficients g. This saves memory but make sure the
!     working mode is: obtain fgco() - use it - obtain again - ...
!
!----------------------------------------------------------------------- 
*                                                                      *
*   This routine evaluates a coefficient                               *
*                                                                      *
*                                K                                     *
*                               g   (IA,IB)                            *
*                                IR                                    *
*                                                                      *
*                                                                      *
*   Here  K  is the multipolarity, IR  is the sequence number of the   *
*   configuration, and  IA and IB are orbital sequence  numbers. See   *
*   I P Grant,  B J McKenzie,  P H Norrington,  D F  Mayers, and N C   *
*   Pyper, Computer Phys Commun 21 (1980) 207-231, Eq (7).             *
*                                                                      *
*   Call(s) to: [LIB92]: CLRX.                                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
! . Two calls to IQ are physically inlined here to reduce the overhead.
! . In the inlined text, the checks for the integer range (ref. 
!    lib92/iq.f for details) are removed since GCO is called only by
!    SETCOF and SETHAM where the input arguments always fall in the
!    suitable range.
!XHH 1997.03.05 
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL FULLA,FULLB,LDBPA
      CHARACTER*2 NH
*
      COMMON/DEBUGA/LDBPA(5)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
!-----------------------------------------------------------------------
!XHH
! New section
! The same COMMON (/fg2cof/) is used to store the coefficients
!
CGG      PARAMETER (NNNWP = 30)
      INTEGER*4 iIQA
      POINTER (PNTRIQ,iIQA(NNNWP,1))
      COMMON/ORB2/NCF,NW,PNTRIQ
      POINTER (pnfgco,fgco(1))
      COMMON/fg2cof/pnfgco

      fac = CLRX (NAK(IA),K,NAK(IB))
      fac2= fac*fac

      DO ir = 1, ncf
         iqa = IBITS(iiqa((ia-1)/4+1,ir),8*MOD(ia-1,4),8)
         iqb = IBITS(iiqa((ib-1)/4+1,ir),8*MOD(ib-1,4),8)
         IF (iqa .EQ. nkj(ia)+1 .OR. iqb .EQ. nkj(ib)+1) THEN
            fgco(ir) = -DBLE(IQA*IQB)*fac2
         ELSE
            fgco(ir) = 0.0D 00
         ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Original section

!      IQA = IQ (IA,IR)
!      IQB = IQ (IB,IR)
!*
!      FULLA = IQA .EQ. NKJ(IA)+1
!      FULLB = IQB .EQ. NKJ(IB)+1
!*
!      IF (FULLA .OR. FULLB) THEN
!         QAB = DBLE (IQA*IQB)
!         FAC = CLRX (NAK(IA),K,NAK(IB))
!         GCO = -QAB*FAC*FAC
!      ELSE
!         GCO = 0.0D 00
!      ENDIF
!*
!*   Debug printout - tens of thousands of times !
!*
!      IF (LDBPA(3) .AND. (ABS (GCO) .GT. 0.0D 00))
!     :      WRITE (99,300) K,NP(IA),NH(IA),NP(IB),NH(IB),GCO,IR
!*
!  300 FORMAT (/'  ',1I2
!     :        /' g    (',1I2,1A2,',',1I2,1A2,') = ',1PD21.14,
!     :        /'  ',1I3/)
!-----------------------------------------------------------------------

      RETURN
      END
