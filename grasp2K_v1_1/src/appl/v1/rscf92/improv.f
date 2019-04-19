************************************************************************
*                                                                      *
      SUBROUTINE IMPROV (EOL,J,lsort)
*                                                                      *
*   Improve the orbital J.                                             *
*                                                                      *
*   Call(s) to: [RSCF92]: DACON, DAMPCK, DAMPOR, LAGCON, MATRIX,       *
*                         NEWCO, ORTHOR, ROTATE, SETCOF, SETLAG,       *
*                         SOLVE, XPOT, YPOT.                           *
*               [LIB92]: ORTHSC, QUAD.                                 *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
*                                                                      *
************************************************************************
*
! CALL ORTHSC (1) --> CALL ORTHSC since it is simply not there
! XHH 1997.02.14
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PCDAMP,PNTRIQ
      POINTER (PCDAMP,CDAMPDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL EOL,FAIL,FIRST,ORTHST,lsort
      CHARACTER*2 NH
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORTHCT/ORTHST
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      PARAMETER (P2    = 2.0D-01,
     :           P005  = 5.0D-03,
     :           P0001 = 1.0D-04)
*
*   C Froese Fischer's IPR and ED1 parameter
*
      DATA IPR /0/
      DATA ED1 /0.0D 00/
      DATA FIRST /.FALSE./

      LOGICAL lcorre
      COMMON/corre/lcorre(NNNW)
      SAVE /corre/
*
      GAMAJ = GAMA(J)
*
*   C Froese Fischer's parameters IPR, ED1, ED2 are set and
*   used in this routine and in DAMPCK
*
    1 ED2 = E(J)
*
*   Set up the exchange potential and arrays XU, XV as appropriate
*
*   Set coefficients for YPOT, XPOT, DACON
*   Compute direct potential, exchange potential
*   Add in Lagrange-multiplier contribution
*   Add in derivative-terms contribution
*
      CALL SETCOF (EOL,J)
      CALL YPOT (J)
      CALL XPOT (J)
      CALL LAGCON (J)
      CALL DACON
*
*   Calculate deferred corrections
*
      CALL DEFCOR (J)
*
*   Solve the Dirac equation
*
      INV = 0
      CALL SOLVE (J,FAIL,INV,JP,NNP)
*
*   Upon failure issue message; take corrective action if possible
*
      IF (FAIL) THEN
         WRITE (*,300) NP(J),NH(J),METHOD(J)
         IF (METHOD(J) .NE. 2) THEN
            METHOD(J) = 2
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
!            CALL ORTHSC (1)
            CALL ORTHSC
            IF (EOL) THEN
               CALL MATRIX
               CALL NEWCO
            ENDIF
            CALL SETLAG (EOL)
            GO TO 1
         ELSE
            WRITE (*,301)
            CALL TIMER (0)
            STOP
         ENDIF
      ENDIF
*
*   Compute norm of radial function
*
      TA(1) = 0.0D 00
      DO 2 I = 2,MTP0
         TA(I) = (P(I)**2+Q(I)**2)*RP(I)
    2 CONTINUE
      MTP = MTP0
      CALL QUAD (DNORM)

!   Determine self-consistency [multiplied by SQRT(UCF(J))]

      CALL CONSIS (J)
*
*   Normalize
*
      DNFAC = 1.0D 00/SQRT (DNORM)
      P0 = P0*DNFAC
      DO 3 I = 1,MTP0
         P(I) = P(I)*DNFAC
         Q(I) = Q(I)*DNFAC
    3 CONTINUE

!   Determine self-consistency [multiplied by SQRT(UCF(J))]
!      CALL CONSIS (J)

*
*   Check if different method should be used or if improvement
*   count should be reduced
*
      DEL1 = ABS (1.0D 00-ED2/E(J))
      IF (METHOD(J) .EQ. 1) THEN
         DEL2 = MAX (ABS (1.0D 00-SQRT (DNORM)),
     :               ABS (DNFAC-1.0D 00))
         IF ((DEL1 .LT. P005) .AND. (DEL2 .GT. P2)) THEN
            METHOD(J) = 2
            GOTO 1
         ENDIF
      ELSE
         IF ((DEL1 .LT. P0001) .AND. (NSIC .GT. 1)) NSIC = NSIC-1
      ENDIF
*
*   Damp the orbital --- if not converged
*
       IF (scnsty(J) .GT. ACCY) THEN
          CALL DAMPCK (IPR,J,ED1,ED2)
          odampj = ABS( odamp(j) )
       ELSE
!          take the whole new orbital
          odampj = 0.D0
       ENDIF
      CALL DAMPOR (J,INV,odampj)

!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.

      IF (ORTHST) THEN
         nwww = nw
         CALL orthy(nwww,J,lsort)
      ENDIF
*
*   Determine self-consistency after damping and orthogonalization
*
!      CALL CONSIS (J)
*
*   Print details of iteration
*
      WRITE (*,302) NP(J),NH(J),E(J),PZ(J),DNORM,METHOD(J),
     :              SCNSTY(J),ODAMPJ,JP,MF(J),INV,NNP
*
      RETURN
*
  300 FORMAT (/' Failure; equation for orbital ',1I2,1A2,
     :         ' could not be solved using method ',1I1)
  301 FORMAT (//' ****** Error in SUBROUTINE IMPROV ******'
     :          /' Convergence not obtained'/)
  302 FORMAT (3X,1I2,1A2,2X,1P,3(1D16.9,2X),2X,1I2,4X,1D11.4,3X,
     :        0P,F6.4,2X,1I4,2X,1I4,3X,1I2,4X,1I2)
*
      END
