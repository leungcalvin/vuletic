************************************************************************
*                                                                      *
      SUBROUTINE NEWCO
*                                                                      *
*   This routine computes the level weights, the generalized occupa-   *
*   tion numbers,  and average energy for  (E)OL calculations;  this   *
*   information and the eigenvectors are then printed out.             *
*                                                                      *
*   Call(s) to: [LIB92]: IQ.                                           *
*               [RSCF92]: CSFWGT, DSUBRS.                              *
*                                         Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
! Another output channel (istde) added to see the energy only
! XHH 
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PIASPA
      POINTER (PIASPA,IASPADUMMY)
      LOGICAL EOL,LDBPG
*
      POINTER (PNTRWT,WT(1))
      POINTER (PWEIGH,WEIGHT(1))
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/DEBUGG/LDBPG(5)
     :      /DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /SCF1/UCF(NNNW)
     :      /SYMA/PIATJP,PIASPA
      COMMON/iounit/istdi,istdo,istde
*
*   Compute weighting factors
*
      SUM = 0.0D 00
      DO 1 J = 1,NCMIN
         WEITJ = WEIGHT(J)
         IF (WEITJ .EQ. -2.0D 00) THEN
            WT(J) = 1.0D 00
         ELSEIF (WEITJ .EQ. -1.0D 00) THEN
            WT(J) = DBLE (IATJPO(J))
         ELSE
            WT(J) = WEIGHT(J)
         ENDIF
         SUM = SUM+WT(J)
    1 CONTINUE
      DO 2 J = 1,NCMIN
         WT(J) = WT(J)/SUM
    2 CONTINUE
*
*   Compute generalised occupation numbers
*
      EOL = .TRUE.
*
      DO 4 J = 1,NW
         SUM = 0.0D 00
         DO 3 I = 1,NCF
            SUM = SUM+DSUBRS (EOL,I,I)*DBLE (IQ (J,I))
    3    CONTINUE
         UCF(J) = SUM
    4 CONTINUE
*
*   Write out level energies and weights
*
      WRITE (*,300)
      SUM = 0.0D 00
      DO 5 J = 1,NCMIN
         II = ICCMIN(J)
         EE = EAV+EVAL(J)
         WRITE (*,301) II,EE,WT(J)
!always         WRITE (istde,301) II,EE,WT(J)
         IF (LDBPG(5)) THEN
            WRITE (99,302)
            WRITE (99,303) (EVEC(I+(J-1)*NCF),I = 1,NCF)
         ENDIF
         SUM = SUM+WT(J)*EE
    5 CONTINUE
      CALL CSFWGT (.TRUE.)
*
*   Write out average energy
*
      IF (NCMIN .GT. 1) WRITE (*,304) SUM
*
*   Write out generalized occupation numbers
*
      WRITE (*,305)
      WRITE (*,303) (UCF(I),I = 1,NW)
*
      RETURN
*
  300 FORMAT (/'Optimise on the following level(s):'/)
  301 FORMAT ('Level ',1I2,4X,'Energy = ',1P,1D19.12,
     :                      4X,'Weight = ',   1D19.12)
  302 FORMAT (/'Configuration mixing coefficients:')
  303 FORMAT (1X,1P,8D15.7)
  304 FORMAT (/'Weighted average energy of these levels = ',1PD18.10)
  305 FORMAT (/'Generalised occupation numbers'/)
*
      END
