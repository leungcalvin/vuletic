************************************************************************
*                                                                      *
      SUBROUTINE PRINTA (ASFA,ASFB,I,J,OMEGA,FACTOR,LINES)
*                                                                      *
*   This  routine  prints the basic oscillator strength  information   *
*   for transitions between level I and level J.                       *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL LTC
      CHARACTER*4 JLBL,LABJ,LABP,JLABI,JLABJ,JPARI,JPARJ
*
      POINTER (PIATJP,IATJPO(1)),(PIASPA,IASPAR(1))
      POINTER(PIATJPII,IATJPOII(1)),(PIASPAII,IASPARII(1))
      POINTER(PIATJPFF,IATJPOFF(1)),(PIASPAFF,IASPARFF(1))
      POINTER (PNTOTB,TOTB(1)),(PNTOTC,TOTC(1))
c
cbieron  numbering of levels in OSCL
c
       POINTER (PNIVEC,IVEC(1))
       POINTER (PNIVECII,IVECII(1))
       POINTER (PNIVECFF,IVECFF(1))
*
      COMMON/CUTO/CUTOFF
     :      /DEF2/C
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMPS,FASI,FBSI
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /SYMA/PIATJP,PIASPA
     :      /OSC2/LK,KK
     :      /OSC4/PNTOTC,PNTOTB
     :      /OSC7/LTC(10)

      COMMON/SYMAII/PIATJPII,PIASPAII
     :      /PRNTII/NVECII,PNIVECII,NVECMXII

      COMMON/SYMAFF/PIATJPFF,PIASPAFF
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF
c
cbieron  numbering of levels in OSCL
c
     :      /PRNT/NVEC,PNIVEC,NVECMX
*
*   Evaluate statistical factors and constants
*
      STFAC = IATJPOII(I)
      DLL1 = DBLE (LK+LK+1)
      STFAC = STFAC/DLL1
*
      OMC = OMEGA/CVAC
c     OMC = OMEGA/C
      FAAU = 2.0D 00*OMC*STFAC/DBLE (IATJPOFF(J))
      FOSC = CVAC*STFAC/OMC
c     FOSC = C*STFAC/OMC
      FBAU = PI*FOSC/OMEGA
      ENG = OMEGA*FACTOR
      IF (LTC(1)) ENG = FACTOR/OMEGA
*
*   J/pi labels for levels
*
      JLABI = LABJ(IATJPOII(I))
      JLABJ = LABJ(IATJPOFF(J))
      IPAR = (IASPARII(I)+3)/2
      JPAR = (IASPARFF(J)+3)/2
      JPARI = LABP(IPAR)
      JPARJ = LABP(JPAR)
*
*   Calculate Einstein A and B coefficients and oscillator strengths
*
      IF (KK .EQ. 0) THEN
*
*   Electric multipoles
*
*   In atomic units
*
         IF (ASFA.LT.0) THEN
            ISIGNA = -1
         ELSE
            ISIGNA = 1
         ENDIF
         IF (ASFB.LT.0) THEN
            ISIGNB = -1
         ELSE
            ISIGNB = 1
         ENDIF

         ACSQ = ASFA**2
         ABSQ = ASFB**2
         AC = ACSQ*FAAU
         AB = ABSQ*FAAU
         BC = ACSQ*FBAU
         BB = ABSQ*FBAU
         OSCC = ACSQ*FOSC
         OSCB = ABSQ*FOSC
         SA=OSCC*1.5/ABS(OMEGA)
         SB=OSCB*1.5/ABS(OMEGA)
*
*   Convert to SI units if option 5 not set
*
         IF (.NOT. LTC(7)) THEN
            AC = AC*FASI
            AB = AB*FASI
            BC = BC*FBSI
            BB = BB*FBSI
         ENDIF
*
*   Accumulate total of A coefficients
*
         TOTC(J) = TOTC(J)+AC
         TOTB(J) = TOTB(J)+AB
*
*   Print information if both AC and AB are greater than CUTOFF
*
         IF ((ABS (AC) .GE. CUTOFF) .AND. (ABS (AB) .GE. CUTOFF))
     :      THEN
c
cbieron  numbering of levels in OSCL
c
c            WRITE (24,300) J,JLABJ,JPARJ,I,JLABI,JPARI,ENG,
!xhh            WRITE (24,300) IVECFF(J),JLABJ,JPARJ,I,JLABI,JPARI,ENG,
            WRITE (24,300) IVECFF(J),JLABJ,JPARJ,IVECII(I),JLABI,JPARI,
     :                     ENG,AC,OSCC,SA,ASFA
C            WRITE (24,*) 'Relative sign',ISIGNA
            WRITE (24,301) AB,OSCB,SB,ASFB
C            WRITE (24,*) 'Relative sign',ISIGNB
            LINES = LINES+3
         ENDIF
*
      ELSE
*
*   Magnetic multipoles
*
*   In atomic units
*
         IF (ASFA.LT.0) THEN
            ISIGNA = -1
         ELSE
            ISIGNA = 1
         ENDIF
         AMS = ASFA**2
         AM = AMS*FAAU
         BM = AMS*FBAU
         OSCM = AMS*FOSC
         SA=OSCM*1.5/ABS(OMEGA)
*
*   Convert to SI units if option 5 not set
*
         IF (.NOT. LTC(7)) THEN
            AM = AM*FASI
            BM = AM*FBSI
         ENDIF
*
*   Accumulate total of A coefficients
*
         TOTC(J) = TOTC(J)+AM
         TOTB(J) = TOTB(J)+AM
*
*   Print information if AM is greater than CUTOFF
*
         IF (ABS (AM) .GE. CUTOFF) THEN
c
cbieron  numbering of levels in OSCL
c
c            WRITE (24,302) J,JLABJ,JPARJ,I,JLABI,JPARI,ENG,
            WRITE (24,302) IVECFF(J),JLABJ,JPARJ,I,JLABI,JPARI,ENG,
     :                      AM,OSCM,SA,ASFA
C            WRITE (24,*) 'Relative sign',ISIGNA
            LINES = LINES+2
         ENDIF
*
      ENDIF
*
      RETURN
*
  300 FORMAT(I3,2A4,I3,2A4,1P,1D15.6,' C',1P,4D15.6)
  301 FORMAT(37X,' B',1P,4D15.6)
  302 FORMAT(I3,2A4,I3,2A4,1P,1D15.6,' M',1P,4D15.6)
*
      END
