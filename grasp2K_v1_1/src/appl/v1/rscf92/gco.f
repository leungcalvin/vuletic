************************************************************************
*                                                                      *
      FUNCTION GCO (K,IR,IA,IB)
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

CGG      PARAMETER (NNNWP = 30)
      INTEGER*4 iIQA
      POINTER (PNTRIQ,iIQA(NNNWP,1))
      COMMON/ORB2/NCF,NW,PNTRIQ

      iqa = IBITS(iiqa((ia-1)/4+1,ir),8*MOD(ia-1,4),8)
      iqb = IBITS(iiqa((ib-1)/4+1,ir),8*MOD(ib-1,4),8)
      IF (iqa .EQ. nkj(ia)+1 .OR. iqb .EQ. nkj(ib)+1) THEN
         FAC = CLRX (NAK(IA),K,NAK(IB))
         GCO = -DBLE(IQA*IQB)*FAC*FAC
      ELSE
         GCO = 0.0D 00
      ENDIF

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
