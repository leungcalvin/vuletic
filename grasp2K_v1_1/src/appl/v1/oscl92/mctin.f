************************************************************************
*                                                                      *
      SUBROUTINE MCTIN (IOPAR,JKP)
*                                                                      *
*   This routine loads  coefficients with parity and  rank specified   *
*   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
*   (+/- 1) and is determined from the sign of  KP(JKP).               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*               [OSCL92]: ALCNSA, ALCNTA, TRSORT.                      *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTJJA,PNTJJB,PNTHB1,PNTHB2,PNTHC1,PNTHC2,PNTHM1,PNTHM2,
Cww     :        PNTRIQ
      POINTER (PNTJJA,JJADUMMY)
      POINTER (PNTJJB,JJBDUMMY)
      POINTER (PNTHB1,HB1DUMMY)
      POINTER (PNTHB2,HB2DUMMY)
      POINTER (PNTHC1,HC1DUMMY)
      POINTER (PNTHC2,HC2DUMMY)
      POINTER (PNTHM1,HM1DUMMY)
      POINTER (PNTHM2,HM2DUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPA,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      CHARACTER*2 NH
*
      DIMENSION TSHELL(NNNW)
*
      POINTER (PLABEL,LABEL(1))
      POINTER (PCOEFF,COEFF(1))
*
      POINTER (PISLDR,ISLDR(1))
      POINTER (PXSLDR,XSLDR(1))
      POINTER (PNTLAB,LAB(1))
      POINTER (PNNPTR,NPTR(1))
      POINTER (PNTRKP,KP(1))
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /FOPARM/ICCUT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /OSC1/PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :            PNTHC1,PNTHC2,PNTHM1,PNTHM2,NSDIM
     :      /OSC2/LK,KK
     :      /OSC3/PXSLDR,PISLDR,NTDIM
     :      /OSC5/NINT,PNTLAB,PNNPTR,NINTEG
     :      /OSC6/NKP,PNTRKP
*
      PARAMETER (CUTOFF = 1.0D-10)
*
CGG      PARAMETER (KEYORB = 121)
*
* unit 99 is connected with the ".dbg" file
* xhh 98-06-29
      PARAMETER (NFILE = 88)
*
*   Open the scratch file to store the MCT coefficients; position file
*   to beginning
*
      print *, '**************** fine before ****************'
      OPEN (NFILE,STATUS = 'unknown', FORM = 'UNFORMATTED')
      print *, '**************** fine after ****************'
      REWIND (NFILE)
*
*   Allocate storage to buffer arrays
*
      NLABEL = 1
      CALL ALLOC (PLABEL,NLABEL,4)
      CALL ALLOC (PCOEFF,NLABEL,8)
*
*   Generate MCT coefficients for given rank/parity combination and
*   store them by CSF on NFILE
*
      LK = ABS (KP(JKP))
      IOPAR = ISIGN (1,KP(JKP))
*
      NMCT = 0
* 
      Print *, '.. Entering MCTIN with rank and Parity :', Lk, IOPAR
      DO 3 IC = 1,NCF
*
         IR = IC
	 IF (MOD (IC,1000) .eq. 0) print *, '  IC =',ic
*
    1    NCR = 0
!	 IF (MOD (IR,1000) .eq. 0) print *, '    IR =',ir
         CALL TNSRJJ (LK,IOPAR,IC,IR,IA,IB,TSHELL)
         IF (IA .NE. 0) THEN
            IF (IA .EQ. IB) THEN
               DO 2 IA = 1,NW
                  IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                     NCR = NCR+1
                     IF (NCR .GT. NLABEL) THEN
                        NEWSIZ = 2*NLABEL
                        CALL RALLOC (PLABEL,NLABEL,NEWSIZ,4)
                        CALL RALLOC (PCOEFF,NLABEL,NEWSIZ,8)
                        NLABEL = NEWSIZ
                     ENDIF
                     LABEL(NCR) = IA*KEYORB+IA
                     COEFF(NCR) = TSHELL(IA)
                  ENDIF
    2          CONTINUE
            ELSE
               IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                  NCR = NCR+1
                  IF (NCR .GT. NLABEL) THEN
                     NEWSIZ = 2*NLABEL
                     CALL RALLOC (PLABEL,NLABEL,NEWSIZ,4)
                     CALL RALLOC (PCOEFF,NLABEL,NEWSIZ,8)
                     NLABEL = NEWSIZ
                  ENDIF
                  LABEL(NCR) = IA*KEYORB+IB
                  COEFF(NCR) = TSHELL(1)
               ENDIF
            ENDIF
         ENDIF
         IF (NCR .GT. 0) THEN
            WRITE (NFILE) IC,IR,NCR
            WRITE (NFILE) (LABEL(I),COEFF(I),I = 1,NCR)
            NMCT = NMCT+NCR
         ENDIF
*
         IF (IR .LT. NCF) THEN
            IF (LFORDR .AND. (IC .GT. ICCUT)) GOTO 3
            IR = IR+1
            GOTO 1
         ENDIF
*
    3 CONTINUE
*
*   Deallocate storage for buffer arrays
*
      CALL DALLOC (PLABEL)
      CALL DALLOC (PCOEFF)
*
      WRITE (*,301) NMCT,LK,IOPAR
*
*   Sort the MCT coefficients by integral labels
*
      CALL TRSORT (NFILE,LDBPA(2))
*
*   Read the data back as required by OSCL conventions
*
      REWIND (NFILE)
*
      M = 0
      K = 0
*
      READ (NFILE) NINT
*
      DO 5 I = 1,NINT
*
         READ (NFILE) LABL,NCSF
*
         M = M+1
         IF (M .GT. NSDIM)
     :      CALL ALCNSA (PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :                   PNTHC1,PNTHC2,PNTHM1,PNTHM2,
     :                   PNTLAB,PNNPTR,NSDIM,2)
         LAB(M) = LABL
*
*   Read configuration pairs and coefficients for this integral
*
    4    IF (NCSF+K .GT. NTDIM) THEN
            CALL ALCNTA (PISLDR,PXSLDR,NTDIM,2)
            GOTO 4
         ENDIF
         NPTR(M) = K
         READ (NFILE) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
         K = K+NCSF
*
    5 CONTINUE
*
*   Close (and hence release) the scratch file
*
      CLOSE (unit=NFILE,status="DELETE")
*
      NPTR(M+1) = K
      NINTEG = M
*
      RETURN
*
  301 FORMAT (///1X,I8,' MCT coefficients generated for rank ',I2,
     :        ' and parity ',I2//)
*
      END

