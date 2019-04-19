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
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
      POINTER (PNTOTB,TOTB(1))
      POINTER (PNTOTC,TOTC(1))
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
*
*   Evaluate statistical factors and constants
*
      STFAC = IATJPO(I)
      DLL1 = DBLE (LK+LK+1)
      STFAC = STFAC/DLL1
*
      OMC = OMEGA/C
      FAAU = 2.0D 00*OMC*STFAC/DBLE (IATJPO(J))
      FOSC = C*STFAC/OMC
      FBAU = PI*FOSC/OMEGA
      ENG = OMEGA*FACTOR
      IF (LTC(1)) ENG = FACTOR/OMEGA
*
*   J/pi labels for levels
*
      JLABI = LABJ(IATJPO(I))
      JLABJ = LABJ(IATJPO(J))
      IPAR = (IASPAR(I)+3)/2
      JPAR = (IASPAR(J)+3)/2
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
         ACSQ = ASFA**2
         ABSQ = ASFB**2
         AC = ACSQ*FAAU
         AB = ABSQ*FAAU
         BC = ACSQ*FBAU
         BB = ABSQ*FBAU
         OSCC = ACSQ*FOSC
         OSCB = ABSQ*FOSC
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
            WRITE (24,300) J,JLABJ,JPARJ,I,JLABI,JPARI,ENG,
     :                      AC,BC,OSCC
            WRITE (24,301) AB,BB,OSCB
            LINES = LINES+3
         ENDIF
*
      ELSE
*
*   Magnetic multipoles
*
*   In atomic units
*
         AMS = ASFA**2
         AM = AMS*FAAU
         BM = AMS*FBAU
         OSCM = AMS*FOSC
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
            WRITE (24,302) J,JLABJ,JPARJ,I,JLABI,JPARI,ENG,
     :                      AM,BM,OSCM
            LINES = LINES+2
         ENDIF
*
      ENDIF
*
      RETURN
*
  300 FORMAT(1X,I4,3X,2A4,I8,3X,2A4,6X,' Coulomb',1P,4D19.7)
  301 FORMAT(41X,'Babushkin',18X,1P,3D19.7/)
  302 FORMAT(1X,I4,3X,2A4,I8,3X,2A4,6X,'Magnetic',1P,4D19.7/)
*
      END
