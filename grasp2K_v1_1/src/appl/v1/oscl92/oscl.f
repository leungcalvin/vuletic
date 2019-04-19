************************************************************************
*                                                                      *
      SUBROUTINE OSCL
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of  data for transitions between multiconfiguration   *
*   Dirac-Fock energy levels.                                          *
*                                                                      *
*   Call(s) to: [OSCL92}: ALCNSA, ALCNTA, BESSJ, CSFM.                 *
*               [LIB92]: ALLOC, DALLOC, MCTIN, PRINTA.                 *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTJJA,PNTJJB,PNTHB1,PNTHB2,PNTHC1,PNTHC2,PNTHM1,PNTHM2,
Cww     :        PNTLAB,PNTRIQ,PNNPTR,PISLDR,PXSLDR,PNTRKP,
Cww     :        PNTRPF,PNTRQF
      POINTER (PNTJJA,JJADUMMY)
      POINTER (PNTJJB,JJBDUMMY)
      POINTER (PNTHB1,HB1DUMMY)
      POINTER (PNTHB2,HB2DUMMY)
      POINTER (PNTHC1,HC1DUMMY)
      POINTER (PNTHC2,HC2DUMMY)                                         
      POINTER (PNTHM1,HM1DUMMY)                                         
      POINTER (PNTHM2,HM2DUMMY)
      POINTER (PNTLAB,LABDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNNPTR,NPTRDUMMY)
      POINTER (PISLDR,ISLDRDUMMY)
      POINTER (PXSLDR,XSLDRDUMMY)
      POINTER (PNTRKP,RKPDUMMY)
      POINTER (PNTRPF,RPFDUMMY)
      POINTER (PNTRQF,RQFDUMMY)

      LOGICAL LTC
      CHARACTER*4 IAU,IEV,ICM,IHZ,IANG,IUNITS
      CHARACTER*2 NH
*
      POINTER (PNTRET,ET(1))
      POINTER (PNTET1,ET1(1))
      POINTER (PNTIPR,IPR(1))
      POINTER (PNIPR1,IPR1(1))
      POINTER (PNNEXT,NEXT(1))
*
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PNTOTB,TOTB(1))
      POINTER (PNTOTC,TOTC(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
*
      COMMON/DEF2/C
     :      /DEF10/AUCM,AUEV,CCMPS,FASI,FBSI
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /OSC1/PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :            PNTHC1,PNTHC2,PNTHM1,PNTHM2,NSDIM
     :      /OSC2/LK,KK
     :      /OSC3/PXSLDR,PISLDR,NTDIM
     :      /OSC4/PNTOTC,PNTOTB
     :      /OSC5/NINT,PNTLAB,PNNPTR,NINTEG
     :      /OSC6/NKP,PNTRKP
     :      /OSC7/LTC(10)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /TITL/IHED,ITIME,IDATE
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      PARAMETER (IAU  = 'Hart',
     :           IEV  = ' eV ',
     :           ICM  = 'Kays',
     :           IHZ  = ' Hz ',
     :           IANG = ' A  ')
      PARAMETER (NCA = 65536)
*
*   Allocate storage
*
      CALL ALLOC (PNTOTB,NVEC,8)
      CALL ALLOC (PNTOTC,NVEC,8)
*
      NVECPR = (NVEC*(NVEC-1))/2
      CALL ALLOC (PNTRET,NVECPR,8)
      CALL ALLOC (PNTET1,NVECPR,8)
      CALL ALLOC (PNTIPR,NVECPR,4)
      CALL ALLOC (PNIPR1,NVECPR,4)
      CALL ALLOC (PNNEXT,NVECPR,4)
*
      CALL ALCNSA (PNTJJA,PNTJJB,PNTHB1,PNTHB2,PNTHC1,PNTHC2,
     :             PNTHM1,PNTHM2,PNTLAB,PNNPTR,NSDIM,1)
      CALL ALCNTA (PISLDR,PXSLDR,NTDIM,1)
*
*   Set up units for printing transition energy
*
      IF     (LTC(1)) THEN
*
*   Print transition energies in Angstroms
*
         FACTOR = AUCM
         FACTOR = 1.0D 08/FACTOR
         IUNITS = IANG
*
      ELSEIF (LTC(2)) THEN
*
*   Print energies in eV
*
         FACTOR = AUEV
         IUNITS = IEV
*
      ELSEIF (LTC(3)) THEN
*
*   Print transition energies in Hartree Atomic Units
*
         FACTOR = 1.0D 00
         IUNITS = IAU
*
      ELSEIF (LTC(4)) THEN
*
*   Print transition energies in Hz
*
         FACTOR = AUCM
         FACTOR = FACTOR*CCMPS
         IUNITS = IHZ
*
      ELSEIF (LTC(5)) THEN
*
*   Print transition energies in Kaysers
*
         FACTOR = AUCM
         IUNITS = ICM
*
      ENDIF
*
*   Initialization for total decay rate
*
      DO 4 I = 1,NVEC
         II = IVEC(I)
         TOTC(II) = 0.0D 00
         TOTB(II) = 0.0D 00
    4 CONTINUE
*
*   Select type of transition
*
*   KK  =  0  for electric multipole.
*       =  1  for magnetic multipole.
*
*   CF:  IOPAR  =  (-1)**N      Electric N-pole
*               =  (-1)**(N+1)  Magnetic N-pole.
*                               N > 0
*
      KKA = 1
      DO 37 JKP = 1,NKP
*
         CALL MCTIN (IOPAR,JKP)
*
         IF (LK .GT. 0) THEN
            IELEC = (-1)**LK
            IF (IELEC .EQ. IOPAR) THEN
               KK = 0
               KKA = 0
            ELSE
               KK = 1
               IELEC = -IELEC
            ENDIF
*
*   Set up list of levels for calculation of oscillator strengths
*   sort list into increasing order of energy if option 6 set
*
            M = 0
            IFIRST = 1
            NEXT(1) = 0
            NVEC1 = NVEC-1
            DO 29 LEVI = 1,NVEC1
               LEVI1 = LEVI+1
               DO 28 LEVF = LEVI1,NVEC
                  LEV1 = IVEC(LEVI)
                  LEV2 = IVEC(LEVF)
                  IF (LEV2 .GE. LEV1) GOTO 24
                  LEV1 = LEV2
                  LEV2 = IVEC(LEVI)
*
*   Check for consistent parity and J
*
   24             ITKPO = LK+LK+1
                  IF (ITRIG (IATJPO(LEV1),IATJPO(LEV2),ITKPO)
     :                .EQ. 0) GOTO 28
                  ITEST = IASPAR(LEV1)*IASPAR(LEV2)*IELEC
                  IF (ITEST .LT. 0) GOTO 28
                  M = M+1
                  IF (M .GT. NVECPR) THEN
                     PRINT *, 'OSCL: Dynamic allocation incorrectly'
                     PRINT *, ' computed: Bug.'
                     STOP
                  ENDIF
                  ET(M) = EVAL(LEV2)-EVAL(LEV1)
                  IPR(M) = LEV1+NCA*LEV2
*
*   Sort level pairs by energy
*
                  IF (LTC(6)) THEN
                     ET1(M) = ET(M)
                     IPR1(M) = IPR(M)
                     L = 0
                     J = IFIRST
   25                IF (ET(M) .GE. ET1(J)) GOTO 26
                     IF (L .NE. 0) GOTO 27
                     NEXT(M) = IFIRST
                     IFIRST = M
                     GOTO 28
   26                L = J
                     J = NEXT(L)
                     IF (J .NE. 0) GOTO 25
   27                NEXT(L) = M
                     NEXT(M) = J
                  ENDIF
*
   28          CONTINUE
   29       CONTINUE
            NTRANS = M
*
*   Unpack sorted list
*
            IF (LTC(6)) THEN
               L = IFIRST
               DO 30 I = 1,NTRANS
                  ET(I) = ET1(L)
                  IPR(I) = IPR1(L)
                  L = NEXT(L)
   30          CONTINUE
            ENDIF
*
*   Calculate and print transition probability data
*
            NLP = 70-8
            LINES = NLP
*
            DO 36 I = 1,NTRANS
*
               IF (LINES .GE. NLP) THEN
                  IF (KK .EQ. 0) THEN
                     WRITE (24,308) LK
                  ELSE
                     WRITE (24,309) LK
                  ENDIF
                  WRITE (24,310)
                  IF (LTC(1)) THEN
                     WRITE (24,311)
                     IF (.NOT. LTC(7)) THEN
                        WRITE (24,312)
                     ELSE
                        WRITE (24,313)
                     ENDIF
                  ELSE
                     WRITE (24,314)
                     IF (.NOT. LTC(7)) THEN
                        WRITE (24,315) IUNITS
                     ELSE
                        WRITE (24,316) IUNITS
                     ENDIF
                  ENDIF
                  LINES = 0
               ENDIF
*
               LEV1 = MOD(IPR(I),NCA)
               LEV2 = IPR(I)/NCA
               OMEGA = -ET(I)
               ARGU = OMEGA/C
               CALL BESSJ (ARGU)
*
*  Calculate oscillator strength between the ASFs
*
               CALL CSFM (ASFA,ASFB,LEV1,LEV2)
               CALL PRINTA (ASFA,ASFB,LEV1,LEV2,OMEGA,FACTOR,
     :                      LINES)
   36       CONTINUE
            WRITE (24,317)
         ENDIF
   37 CONTINUE
*
*   Print lifetimes and widths of levels
*
      NLP = 70-10
      LINES = NLP
      DO 42 I = 1,NVEC
         IF (LINES .GE. NLP) THEN
            WRITE (24,318)
            LINES = 0
         ENDIF
         J = IVEC(I)
         TTC = ABS (TOTC(J))
         TTB = ABS (TOTB(J))
         IF ((TTC .NE. 0.0D 00) .AND.
     :       (TTB .NE. 0.0D 00)) THEN
            IF (.NOT. LTC(7)) THEN
               TCCM = TTC/CCMPS
               TBCM = TTB/CCMPS
               TCAU = TCCM/AUCM
               TBAU = TBCM/AUCM
               TCSEC = 1.0D 00/TTC
               TBSEC = 1.0D 00/TTB
            ELSE
               TCAU = TTC
               TBAU = TTB
               TCCM = TTC*AUCM
               TBCM = TTB*AUCM
               TCSEC = 1.0D 00/(TCCM*CCMPS)
               TBSEC = 1.0D 00/(TBCM*CCMPS)
            ENDIF
            TCEV = TCAU*AUEV
            TBEV = TBAU*AUEV
            IF (KKA .NE. 1) THEN
               WRITE (24,319) J,TCSEC,TCAU,TCCM,TCEV
               WRITE (24,320)   TBSEC,TBAU,TBCM,TBEV
               LINES = LINES+3
            ELSE
               WRITE (24,321) J,TCSEC,TCAU,TCCM,TCEV
               LINES = LINES+2
            ENDIF
         ENDIF
   42 CONTINUE
*
*   Deallocate storage; this is local to OSCL
*
  999 CALL DALLOC (PNTOTB)
      CALL DALLOC (PNTOTC)
*
      CALL DALLOC (PNTRET)
      CALL DALLOC (PNTET1)
      CALL DALLOC (PNTIPR)
      CALL DALLOC (PNIPR1)
      CALL DALLOC (PNNEXT)
*
      CALL DALLOC (PNTRKP)
*
      CALL ALCNSA (PNTJJA,PNTJJB,PNTHB1,PNTHB2,PNTHC1,PNTHC2,
     :             PNTHM1,PNTHM2,PNTLAB,PNNPTR,NSDIM,3)
      CALL ALCNTA (PISLDR,PXSLDR,NTDIM,3)
*
*   This was allocated in LOAD
*
      CALL DALLOC (PNTRPF)
      CALL DALLOC (PNTRQF)
      CALL DALLOC (PNIVEC)
      CALL DALLOC (PIATJP)
      CALL DALLOC (PIASPA)
      CALL DALLOC (PNEVAL)
      CALL DALLOC (PNEVEC)
*
*   Close all files
*
      CLOSE (24)
*
      RETURN
*
  302 FORMAT (/' ***** Warning *****')
  303 FORMAT (//' ***** Error in OSCL *****')
  307 FORMAT (/' Dynamic allocation computed incorrectly: Bug.')
  308 FORMAT (' Electric 2**(',I2,')-pole transitions')
  309 FORMAT (' Magnetic 2**(',I2,')-pole transitions')
  310 FORMAT (1X,33('='))
  311 FORMAT (/'   Upper state        Lower state  ',8X,'Gauge',8X,
     : 'Wavelength',13X,'Einstein coefficients',13X,'Oscillator')
  312 FORMAT(81X,'-1',15X,'3 -2 -1',/' Level  J Parity',4X,'Level  J ',
     : 'Parity',21X,'(Angstroms)',10X,'A (s  )',9X,'gB (m s  J  )',7X,
     : 'strength gf'/)
  313 FORMAT(' Level  J Parity',4X,'Level  J Parity',21X,'(Angstroms)',
     : 10X,'A (au)',13X,'gB (au)',10X,'strength gf'/)
  314 FORMAT (/'   Upper state        Lower state  ',8X,'Gauge',10X,
     : 'Energy',15X,'Einstein coefficients',13X,'Oscillator')
  315 FORMAT(81X,'-1',15X,'3 -2 -1',/' Level  J Parity',4X,'Level  J ',
     : 'Parity',23X,'(',A4,')',13X,'A (s  )',9X,'gB (m s  J  )',7X,
     : 'strength gf'/)
  316 FORMAT(' Level  J Parity',4X,'Level  J Parity',23X,'(',A4,')',13X,
     : 'A (au)',13X,'gB (au)',10X,'strength gf'/)
  317 FORMAT (/1X,124('+'))
  318 FORMAT (' Radiative lifetimes and widths'
     :       /' =============================='
     :      //' Level',6X,'Gauge',12X,'Lifetime',33X,'Width'
     :       /' -----',6X,'-----',12X,'--------',12X,
     :        '---------------------------------------------',
     :   /29X,'seconds',13X,'Hartrees',12X,'Kaysers',16X,'eV'/)
  319 FORMAT (1X,I4,6X,'Coulomb: ',1P,4D20.7)
  320 FORMAT (10X,'Babushkin:',1P,4D20.7/)
  321 FORMAT (1X,I4,5X,'Magnetic: ',1P,4D20.7/)
*
      END
