************************************************************************
*                                                                      *
      SUBROUTINE MCTIN (IOPAR,JKP,NAME)
*                                                                      *
*   This routine loads  coefficients with parity and  rank specified   *
*   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
*   (+/- 1) and is determined from the sign of  KP(JKP).               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*               [OSCL92]: ALCNSA, ALCNTA, TRSORT.                      *
*                                                                      *
*                                         Last revision: 28 Dec 1992   *
*   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (KEYORB = 121)

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
      LOGICAL LDBPA,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,AVAIL
      CHARACTER*2 NH
      CHARACTER*24 NAME(2)
*
      DIMENSION TSHELL(NNNW)
*
      POINTER (PLABEL,LABEL(1)),(PCOEFF,COEFF(1))
*
      POINTER (PISLDR,ISLDR(1)),(PXSLDR,XSLDR(1))
      POINTER (PNTLAB,LAB(1)),(PNNPTR,NPTR(1))
      POINTER (PNTRKP,KP(1))
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /FOPARM/ICCUT
     :      /OFFD/NOFFD1,NOFFD2
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
      PARAMETER (NFILE = 93)
      PARAMETER (NFILE1 = 237)
*
*   Check if angular data is available on file
*
      NFILE2 = NFILE1 +JKP
      CALL ANGDATA(NAME,AVAIL,JKP,NFILE2)
*
*   If angular data is not available open the scratch file to store the MCT
*   coefficients; position file
*   to beginning
*
      IF (.NOT.AVAIL) THEN
        OPEN (NFILE,STATUS = 'unknown', FORM = 'UNFORMATTED')
        REWIND (NFILE)
      ENDIF
*
      LK = ABS (KP(JKP))
      IOPAR = ISIGN (1,KP(JKP))
*
      NMCT = 0
*
*   
*   If angular data is not available
*   Generate MCT coefficients for given rank/parity combination and
*   store them by CSF on NFILE
*
      IF (.NOT.AVAIL) THEN
*
*   Allocate storage to buffer arrays
*
        NLABEL = 32
        CALL ALLOC (PLABEL,NLABEL,4)
        CALL ALLOC (PCOEFF,NLABEL,8)

        DO 3 IC = 1,NCF
*
           IR = IC
*
    1      NCR = 0

*                                                                        
*   In many case one is interested only in M1 and E2 transitions between
*   levels with different J values. If this is the case then the do check
*   on the J quantum numbers of the CSFs before calling TNSRJJ.         
*                                                                       
           IF (KP(JKP).EQ.1.AND.NOFFD1.EQ.1) THEN                         
             IF (ITJPO (IC).EQ.ITJPO (IR)) GOTO 13                        
           ENDIF                                                          
           IF (KP(JKP).EQ.2.AND.NOFFD2.EQ.1) THEN                         
             IF (ITJPO (IC).EQ.ITJPO (IR)) GOTO 13                        
           ENDIF
           if(ispar(ic)*ispar(ir)*iopar.ne.1.
     &        or.itrig(itjpo(ic),itjpo(ir),2*lk+1).ne.1) go to 13
c          if(ichkq1(IC,IR).eq.0) go to 13
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
    2            CONTINUE
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
c           print *, 'IC,IR,IA,IB,COEFF(I)'
c         do i=1,ncr
c           write(*,'(4i3,e14.7)') IC,IR,label(i)/keyorb,
c     :  mod(label(i),keyorb),COEFF(i)
c        enddo
              NMCT = NMCT+NCR
           ENDIF
*
*                                                                       
   13      CONTINUE                                                       
                                                                        
           IF (IR .LT. NCF) THEN
              IF (LFORDR .AND. (IC .GT. ICCUT)) GOTO 3
              IR = IR+1
              GOTO 1
           ENDIF
*
    3   CONTINUE
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
        CALL TRSORT (NAME,NFILE,NFILE2,LDBPA(2),JKP)
*
*   Read the data back as required by OSCL conventions
*
        REWIND (NFILE)
      ENDIF

      M = 0
      K = 0
*
      IF (AVAIL) THEN
        READ (NFILE2) NINT
      ELSE
        READ (NFILE) NINT
      ENDIF
*
      DO 5 I = 1,NINT
*
         IF (AVAIL) THEN
           READ (NFILE2) LABL,NCSF
         ELSE
           READ (NFILE) LABL,NCSF
         ENDIF
*
         M = M+1
cbieron
         IF (M .GE. NSDIM)
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
         IF (AVAIL) THEN
           READ (NFILE2) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
         ELSE
           READ (NFILE) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
*           write(*,*) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
         ENDIF
         K = K+NCSF
*
    5 CONTINUE
*
*   Close (and hence release) the scratch file
*
      IF (AVAIL) THEN
        CLOSE (NFILE2)
      ELSE
        CLOSE (unit=NFILE,status="DELETE")
      ENDIF
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
