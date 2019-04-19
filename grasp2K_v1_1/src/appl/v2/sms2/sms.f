************************************************************************
*                                                                      *
      SUBROUTINE SMS
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the  sms parameter, the electron density at the    *
*   origin.
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, GETYN          *
*                        ITJPO, RKCO, TNSRJJ                           *
*               [SMS92]: RINTISO, RINTDENS, VINTI                      *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
*                                         Last revision: 10 Nov 1995   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ,
Cww     :        PINDTE,PVALTE
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PINDTE,INDTEDUMMY)
      POINTER (PVALTE,VALTEDUMMY)
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      CHARACTER*2 CK,NH
      LOGICAL GETYN,FIRSTT,LDBPA,VSH,NUCDE,SMSSH,YES,AVAIL
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION VINT(NNNW,NNNW),TSHELL(NNNW)
      DIMENSION DINT1(NNNW,NNNW),DINT2(NNNW,NNNW),DINT3(NNNW,NNNW)
      DIMENSION DINT4(NNNW,NNNW),DINT5(NNNW,NNNW),DINT6(NNNW,NNNW)
*
      POINTER (PNSMS,SMSC(1))
      POINTER (PNDENS1,DENS1(1))
      POINTER (PNDENS2,DENS2(1))
      POINTER (PNDENS3,DENS3(1))
      POINTER (PNDENS4,DENS4(1))
      POINTER (PNDENS5,DENS5(1))
      POINTER (PNDENS6,DENS6(1))
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNIVEC,IVEC(1))
      POINTER (PIATJP,IATJPO(1)),(PIASPA,IASPAR(1))
*
      EXTERNAL COR,CORD
      COMMON/DEBUGA/LDBPA(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF3/EMPAM,RBCM
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /DEF11/FMTOAU,B1
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /NPAR/PARM(2),NPARM
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /OPT6/NTC(10)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /TEILST/NDTEA,NTEI,PINDTE,PVALTE,FIRSTT
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /SMS1/PNSMS,PNDENS1,PNDENS2,PNDENS3,PNDENS4,PNDENS5,PNDENS6
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
*   Allocate storage for local arrays
*
      CALL ALLOC (PNSMS,NVEC,8)
      CALL ALLOC (PNDENS1,NVEC,8)
      CALL ALLOC (PNDENS2,NVEC,8)
      CALL ALLOC (PNDENS3,NVEC,8)
      CALL ALLOC (PNDENS4,NVEC,8)
      CALL ALLOC (PNDENS5,NVEC,8)
      CALL ALLOC (PNDENS6,NVEC,8)
*
*   Initialise
*
      DO 1 I = 1,NVEC
         SMSC(I) = 0.0D 00
         DENS1(I) = 0.0D 00
         DENS2(I) = 0.0D 00
         DENS3(I) = 0.0D 00
         DENS4(I) = 0.0D 00
         DENS5(I) = 0.0D 00
         DENS6(I) = 0.0D 00
    1 CONTINUE

      VSH   = .TRUE.
      SMSSH = .TRUE.
*
*   Calculate all integrals needed for the volume shift calc.
*
      IF (VSH) THEN
         DO 5 I = 1,NW
            DO 4 J = 1,NW
               IF (NAK(I).EQ.NAK(J)) THEN
                  DINT1(I,J) = RINTDENS(I,J)
                  DINT2(I,J) = RINTI(I,J,1)
                  DINT3(I,J) = RINT(I,J,1)
                  DINT4(I,J) = RINT(I,J,2)
                  DINT5(I,J) = RINT(I,J,-1)
                  DINT6(I,J) = RINT(I,J,-2)
               ELSE
                  DINT1(I,J) = 0.0D 00
                  DINT2(I,J) = 0.0D 00
                  DINT3(I,J) = 0.0D 00
                  DINT4(I,J) = 0.0D 00
                  DINT5(I,J) = 0.0D 00
                  DINT6(I,J) = 0.0D 00
               ENDIF
    4       CONTINUE
    5    CONTINUE
      ENDIF
*
*   Calculate and save the Vinti integrals
*
      IF (SMSSH) THEN
         DO 7 I = 1,NW
            DO 6 J = 1,NW
               IF (I.NE.J) THEN
                  VINT(I,J) = VINTI(I,J)
               ELSE
                  VINT(I,J) = 0.0D 00  
               ENDIF
    6       CONTINUE
    7    CONTINUE
      ENDIF
*
*   See if the appropriate angular data is available. If so,
*   then read the angular files and perform the calculation.
*   
      IF (SMSSH) THEN
         CALL SETMCP(AVAIL)
         IF (AVAIL) THEN
            CALL SMSMCP(VINT)
         ELSE
            CALL SMSNEW(VINT) 
         ENDIF
      ENDIF
      IF (VSH) THEN
         CALL SETMCP(AVAIL)
         IF (AVAIL) THEN
            CALL DENSMCP(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6)
         ELSE
            CALL DENSNEW(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6)
         ENDIF
      ENDIF
*
*   Printouts
*
CPJ      WRITE (24,301) CUTOFF
      WRITE (24,302)
      DO 14 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  SMSC(I)
   14 CONTINUE
      WRITE (24,307)
      DO 16 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS1(I)
   16 CONTINUE
      WRITE (24,308)
      DO 18 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS2(I)
   18 CONTINUE
      WRITE (24,309)
      DO 20 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS3(I)
   20 CONTINUE
      WRITE (24,310)
      DO 22 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS4(I)
   22 CONTINUE
      WRITE (24,311)
      DO 24 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS5(I)
   24 CONTINUE
      WRITE (24,312)
      DO 26 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
     :                  DENS6(I)
   26 CONTINUE
*
*   Dealloc
*
      CALL DALLOC (PNSMS)
      CALL DALLOC (PNDENS1)
      CALL DALLOC (PNDENS2)
      CALL DALLOC (PNDENS3)
      CALL DALLOC (PNDENS4)
      CALL DALLOC (PNDENS5)
      CALL DALLOC (PNDENS6)
      RETURN
*
  301 FORMAT (//' CUTOFF set to ',1PD22.15)
  302 FORMAT (//' Level  J Parity  Specific mass shift (au) '/)
  303 FORMAT (1X,I3,5X,2A4,3X,D20.10)
  307 FORMAT (//' Electron density in atomic units'
     :        //' Level  J Parity',8X,'DENS (a.u.)'/)
  308 FORMAT (//' Kinetic energy '
     :        //' Level  J Parity',8X,'T (a.u.)'/)
  309 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r> (a.u.)'/)
  310 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r2> (a.u.)'/)
  311 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r-1> (a.u.)'/)
  312 FORMAT (//' Radial expectationvalue'
     :        //' Level  J Parity',8X,'<r-2> (a.u.)'/)

*
      END
