************************************************************************
*                                                                      *
      SUBROUTINE ENDSUM
*                                                                      *
*   Generates the last part of  rscf92.sum  (on stream 24).            *
*                                                                      *
*   Call(s) to: [LIB92]: ENGOUT, RINT.                                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*2 NH
*
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /SCF1/UCF(NNNW)
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /SYMA/PIATJP,PIASPA
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Write out the orbital properties
*
      WRITE (24,301)
      DO 1 I = 1,NW
         WRITE (24,302) NP(I),NH(I),E(I),PZ(I),
     :                  GAMA(I),PF(2,I),QF(2,I),MF(I)
    1 CONTINUE
*
      WRITE (24,303)
      DO 2 I = 1,NW
         WA = RINT (I,I,-1)
         WB = RINT (I,I, 1)
         WC = RINT (I,I, 2)
         WRITE (24,304) NP(I),NH(I),WA,WB,WC,UCF(I),SCNSTY(I)
    2 CONTINUE
*
      IF (NCMIN .NE. 0) THEN
         CALL ENGOUT (EAV,EVAL,IATJPO,IASPAR,ICCMIN,NCMIN,MODE)
         CALL CSFWGT (.FALSE.)
      ENDIF
*
      CLOSE (24)
*
      RETURN
*
  301 FORMAT (/'Radial wavefunction summary:'
     :       //' Subshell',11X,'e',20X,'p0',18X,
     :         'gamma',19X,'P(2)',18X,'Q(2)',10X,'MTP'
     :        /)
  302 FORMAT (3X,1I2,1A2,1X,1P,5(3X,1D19.12),3X,1I3)
  303 FORMAT (/21X,'-1',42X,'2',15X,'Generalised',14X,'Self-'
     :        /' Subshell',7X,'<  r     >',13X,'<  r  >',14X,
     :         '<  r    >',13X,'occupation',11X,'consistency'/)
  304 FORMAT (3X,1I2,1A2,1X,1P,5(3X,1D19.12))
*
      END
