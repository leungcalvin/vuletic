************************************************************************
*                                                                      *
      SUBROUTINE ORBOUT
*                                                                      *
*   Write all subshell radial wavefunctions to the  .rwf  file.        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Position the file
*

! To save the intermediate results, the output file of the radial 
! wavefunctions is openned and closed for each iteration.
! A CLOSE statement is also addded in setrwg.f where the rwfn.out file
! was originally openned.
! XHH 1997.02.17
!      REWIND (23)
!      READ (23)
      OPEN(23, FILE='rwfn.out', STATUS='UNKNOWN',FORM='UNFORMATTED')
      WRITE (23) 'G92RWF'
*
*   Write out the radial wavefunctions
*
      DO 1 J = 1,NW
         MFJ = MF(J)
         WRITE (23) NP(J),NAK(J),E(J),MFJ
         WRITE (23) PZ(J),(PF(I,J),I = 1,MFJ),(QF(I,J),I = 1,MFJ)
         WRITE (23) (R(I),I = 1,MFJ)
    1 CONTINUE
!
      CLOSE(23)
*
      RETURN
      END
