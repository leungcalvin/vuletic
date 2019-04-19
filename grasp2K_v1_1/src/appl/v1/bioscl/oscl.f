************************************************************************
*                                                                      *
      SUBROUTINE OSCL(NAME)
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of  data for transitions between multiconfiguration   *
*   Dirac-Fock energy levels.                                          *
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
      LOGICAL LTC,AVAIL,lsame
      CHARACTER*4 IAU,IEV,ICM,IHZ,IANG,IUNITS
      CHARACTER*2 NH
      CHARACTER*24 NAME(2)
*
      POINTER (PNTRET,ET(1)),(PNTET1,ET1(1)),
     :        (PNTIPR,IPR(1)),(PNIPR1,IPR1(1)),
     :        (PNNEXT,NEXT(1))
*
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PNTOTB,TOTB(1)),(PNTOTC,TOTC(1))
      POINTER (PIATJP,IATJPO(1)),(PIASPA,IASPAR(1))
*
      COMMON/DEF2/C
     :      /DEF9/CVAC,PI
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
      POINTER(PNEVALII,EVALII(1))
      POINTER(PNEVECII,EVECII(1))
      POINTER(PNIVECII,IVECII(1))
      POINTER(PIATJPII,IATJPOII(1)),(PIASPAII,IASPARII(1))

      POINTER(PNEVALFF,EVALFF(1))
      POINTER(PNEVECFF,EVECFF(1))
      POINTER(PNIVECFF,IVECFF(1))
      POINTER(PIATJPFF,IATJPOFF(1)),(PIASPAFF,IASPARFF(1))

      COMMON/DEF1II/EMNII,IONCTYII,NELECII,ZII
     :      /EIGVALII/EAVII,PNEVALII
     :      /EIGVECII/PNEVECII
     :      /ORB2II/NCFII,NWII
     :      /PRNTII/NVECII,PNIVECII,NVECMXII
     :      /SYMAII/PIATJPII,PIASPAII

      COMMON/DEF1FF/EMNFF,IONCTYFF,NELECFF,ZFF
     :      /EIGVALFF/EAVFF,PNEVALFF
     :      /EIGVECFF/PNEVECFF
     :      /ORB2FF/NCFFF,NWFF
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF
     :      /SYMAFF/PIATJPFF,PIASPAFF

      PARAMETER (IAU  = 'Hart',
     :           IEV  = ' eV ',
     :           ICM  = 'Kays',
     :           IHZ  = ' Hz ',
     :           IANG = ' A  ')
      PARAMETER (NCA = 65536)
      lsame = name(1).eq.name(2)
*
*   Allocate storage
*
      CALL ALLOC (PNTOTB,NVECFF,8)
      CALL ALLOC (PNTOTC,NVECFF,8)
*
      NVECPR = NVECII*NVECFF
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
*   Make a connection between the orbitals of the
*   merged list and the initial and final state lists
*
      CALL CONNECT
*
*   Set up units for printing transition energy
*
      IF (LTC(1)) THEN
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
      DO I = 1,NVECFF
        TOTC(I) = 0.0D 00
        TOTB(I) = 0.0D 00
      ENDDO
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

      DO JKP = 1,NKP
*
        CALL MCTIN (IOPAR,JKP,NAME)
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
*
          DO LEVII = 1,NVECII
            DO LEVFF = 1,NVECFF
C              IF ((EVALFF(LEVFF)+EAVFF).LT.(EVALII(LEVII)+EAVII)) THEN
C                WRITE(*,*) 'Initial state higher than final state'
C                STOP
C              ENDIF
*
*   Check for consistent parity and J
*
              ITKPO = LK+LK+1
              IF (ITRIG (IATJPOII(LEVII),IATJPOFF(LEVFF),ITKPO)
     :                .EQ. 0) GOTO 28
              ITEST = IASPARII(LEVII)*IASPARFF(LEVFF)*IELEC
              IF (ITEST .LT. 0) GOTO 28
*
*   Calculate and print transition probability data
*
              NLP = 70-8
              LINES = NLP
*
              IF (LINES .GE. NLP) THEN
                LINES = 0
              ENDIF
*
              M = LEVFF+NVECFF*(LEVII-1)
c             M = LEVFF+NVECII*(LEVII-1)
              ET(M) = (EVALFF(LEVFF)+EAVFF-EVALII(LEVII)-EAVII)
              if(lsame.and.(et(m).le.0.0)) cycle
Cww              IF (ET(M).GT.0) THEN
Cww                NHIGH(LEVII,LEVFF) = 1
Cww              ELSE
Cww                NHIGH(LEVII,LEVFF) = -1
Cww                ET(M) = -ET(M)
Cww              ENDIF
Cww              OMEGA = ET(M)
              OMEGA = -ET(M)
              ARGU = OMEGA/C
              CALL BESSJ (ARGU)
*
*  Calculate oscillator strength between the ASFs
*
              CALL CSFM (ASFA,ASFB,LEVII,LEVFF)
              CALL PRINTA (ASFA,ASFB,LEVII,LEVFF,OMEGA,FACTOR,LINES)
C              WRITE (24,317)
   28       ENDDO
          ENDDO
        ENDIF
      ENDDO
*
*   Print lifetimes and widths of levels
*
C      NLP = 70-10
C      LINES = NLP
C      DO I = 1,NVECFF
C        IF (LINES .GE. NLP) THEN
C          WRITE (24,318)
C          LINES = 0
C        ENDIF
C        J = I
C        TTC = ABS (TOTC(J))
C        TTB = ABS (TOTB(J))
C        IF ((TTC .NE. 0.0D 00) .AND.(TTB .NE. 0.0D 00)) THEN
C          IF (.NOT. LTC(7)) THEN
C            TCCM = TTC/CCMPS
C            TBCM = TTB/CCMPS
C            TCAU = TCCM/AUCM
C            TBAU = TBCM/AUCM
C            TCSEC = 1.0D 00/TTC
C            TBSEC = 1.0D 00/TTB
C          ELSE
C            TCAU = TTC
C            TBAU = TTB
C            TCCM = TTC*AUCM
C            TBCM = TTB*AUCM
C            TCSEC = 1.0D 00/(TCCM*CCMPS)
C            TBSEC = 1.0D 00/(TBCM*CCMPS)
C          ENDIF
C          TCEV = TCAU*AUEV
C          TBEV = TBAU*AUEV
C          IF (KKA .NE. 1) THEN
C            WRITE (24,319) IVECFF(J),TCSEC
C            WRITE (24,320)   TBSEC
C            LINES = LINES+3
C          ELSE
C            WRITE (24,321) J,TCSEC
C            LINES = LINES+2
C          ENDIF
C        ENDIF
C      ENDDO
*
*   Deallocate storage; this is local to OSCL
*
      CALL DALLOC (PNTOTB)
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
      if(PNTRPF.ne.0) CALL DALLOC (PNTRPF)
      if(PNTRQF.ne.0) CALL DALLOC (PNTRQF)
      if(PNIVEC.ne.0) CALL DALLOC (PNIVEC)
      if(PIATJP.ne.0) CALL DALLOC (PIATJP)
      if(PIASPA.ne.0) CALL DALLOC (PIASPA)
      if(PNEVAL.ne.0) CALL DALLOC (PNEVAL)
      if(PNEVEC.ne.0) CALL DALLOC (PNEVEC)
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
  308 FORMAT (//' Electric 2**(',I2,')-pole transitions')
  309 FORMAT (//' Magnetic 2**(',I2,')-pole transitions')
  310 FORMAT (1X,33('='))
  311 FORMAT (/'   Upper state        Lower state  ',8X,'Gauge',8X,
     : 'Wavelength',13X,'Einstein coefficients',13X,'Oscillator')
  312 FORMAT(81X,'-1',15X,'3 -2 -1',/' Level  J Parity',4X,'Level  J ',
     : 'Parity',21X,'(Angstroms)',10X,'A (s  )',9X,'gB (m s  J  )',7X,
     : 'strength gf'/)
  313 FORMAT(' Level  J Parity',4X,'Level  J Parity',21X,'(Angstroms)',
     : 10X,'A (au)',13X,'gB (au)',10X,'strength gf'/)
  314 FORMAT (/' Upper       Lower ')   
  315 FORMAT(' Lev  J P',3X,'Lev  J P',
     : 7X,'E (',A4,')',9X,'A (s-1)',10X,
     : 'gf',12X,'S',14X,'M')
  316 FORMAT(' Level  J Parity',4X,'Level  J Parity',23X,'(',A4,')',13X,
     : 'A (au)',13X,'gB (au)',10X,'strength gf'/)
  317 FORMAT (/1X,124('+'))
  318 FORMAT (//' Radiative lifetimes '
     :       /' ======================='
     :      //' Level      Lifetime  s (-1)')
  319 FORMAT (1X,I4,6X,'Coulomb: ',1P,1D20.7)
  320 FORMAT (10X,'Babushkin:',1P,1D20.7/)
  321 FORMAT (1X,I4,5X,'Magnetic: ',1P,1D20.7/)
*
      END
