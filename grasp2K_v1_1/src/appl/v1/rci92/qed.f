************************************************************************
*                                                                      *
      SUBROUTINE QED (jstate,SLFINT)
*                                                                      *
*   This  routine estimates corrections to the  energy levels due to   *
*   self-energy.                                                       *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SCREEN.      *
*               [RCI92]: FZALF, HOVLAP.                                *
*                                                                      *
*                                           Last update: 30 Oct 1992   *
*   Modified by Xinghong He                 Last update: 24 Jun 1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ,PNIVEC
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNIVEC,IVECDUMMY)
      CHARACTER*2 NH, npchar*1, nakchar
*
      POINTER (PNTUCF,UCF(1))
*
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      DIMENSION PTEMP(NNNP),QTEMP(NNNP),SLFINT(NNNW)
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF9/CVAC,PI
     :      /EIGVEC/PNEVEC
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /HORB/PH(NNNP),QH(NNNP)
     :      /NPAR/PARM(2),NPARM
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge. 
!
      maxiter = 20
*
*   Determine `generalised occupation numbers'
*
      CALL ALLOC (PNTUCF,NW,8)
!
! Modified so that UCFJ describes the current eigenstate
!
      DO 4 J = 1,NW
         UCFJ = 0.0D 00
!         DO 3 I = 1,NVEC
          I = jstate
            DO 2 II = 1,NCF
               UCFJ = UCFJ+DBLE (IQ (J,II))
     :                   *EVEC(II+(I-1)*NCF)**2
    2       CONTINUE
!    3    CONTINUE
         UCF(J) = UCFJ/DBLE (NCF)
    4 CONTINUE
*
      DO 14 J = 1,NW
*
         NPJ = NP(J)
*
         IF (NPJ .LE. 8) THEN
*
*   Only orbitals with principal quantum number 8 or less can
*   be treated by this section of code
*
            KAPPA = NAK(J)
*
*   Begin by transferring the function to a temporary array
*
            MFJ = MF(J)
*
            PTEMP(1) = 0.0D 00
            QTEMP(1) = 0.0D 00
            DO 5 I = 2,MFJ
               PTEMP(I) = PF(I,J)
               QTEMP(I) = QF(I,J)
    5       CONTINUE
*
*   Determine an effective charge by finding the hydrogenic orbital
*   with the same symmetry with the maximum overlap with the
*   present orbital
*
*   Begin by bracketing the charge around the maximum overlap
*
            ZEFF = Z-SCREEN (J,UCF)
!
! Find the approximate zeff
! maxiter is an arbitarily imposed maximum number foriteration
!
            ZTRY = ZEFF
            DO i = 1, maxiter
               OVRLAP = HOVLAP (PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZTRY)
               OVRLAP = ABS (OVRLAP)
               IF (OVRLAP .LT. 0.5D 00) THEN
                  ztry = ztry + zeff
                  IF (ztry .GE. c) THEN
                     zeff = c-zeff
!                     ... At this exit, zeff does not depend on 
!                     ... the iteration ...
                     EXIT
                  ENDIF
               ELSE
                  ZEFF = ZTRY
!                  ... This would be the normal exit ...
                  EXIT
               ENDIF
            ENDDO

!            PRINT *, jstate, np(j),nh(j),nak(j),' first=',i

            ZINC = 0.1D 00*ZEFF
*
*   Increase the charge until an upper limit on the charge is found
*
            ZMAX = ZEFF
            DO i = 1, maxiter
               ZMAX = ZMAX+ZINC
               IF (ZMAX .LT. C) THEN
                  OVRMAX = HOVLAP (PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZMAX)
                  OVRMAX = ABS (OVRMAX)
                  IF (OVRMAX .LE. OVRLAP) EXIT
               ELSE
                  ZMAX = C
                  EXIT
               ENDIF
            ENDDO
!            PRINT *, jstate, np(j),nh(j),nak(j),' second=',i

*
*   Decrease the charge until a lower limit on the charge is found
*
            ZMIN = ZEFF
            DO i = 1, maxiter
               ZMIN = ZMIN-ZINC
               IF (ZMIN .GT. 0.1D 00) THEN
                  OVRMIN = HOVLAP (PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZMIN)
                  OVRMIN = ABS (OVRMIN)
                  IF (OVRMIN .GT. OVRLAP) CYCLE
               ELSE
                  ZMIN = 0.1D 00
                  EXIT
               ENDIF
            ENDDO
!            PRINT *, jstate, np(j),nh(j),nak(j),' third=',i

*
*   Now find the maximum
*
            DO i = 1, maxiter
               ZMXTRY = 0.5D 00*(ZMAX+ZEFF)
               OVRMXT = HOVLAP (PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZMXTRY)
               OVRMXT = ABS (OVRMXT)
               ZMNTRY = 0.5D 00*(ZMIN+ZEFF)
               OVRMNT = HOVLAP (PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZMNTRY)
               OVRMNT = ABS (OVRMNT)
               IF (OVRMXT .LT. OVRLAP) THEN
                  ZMAX = ZMXTRY
               ENDIF
               IF (OVRMNT .LT. OVRLAP) THEN
                  ZMIN = ZMNTRY
               ENDIF
               ZEFF = 0.5D 00*(ZMAX+ZMIN)
               OVRLAP = HOVLAP (PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZEFF)
               IF (ABS (ZMAX-ZMIN)/ZEFF .GT. 1.0D-02) THEN
                  OVRLAP = ABS (OVRLAP)
               ELSE
                  EXIT
               ENDIF
            ENDDO
!            PRINT *, jstate, np(j),nh(j),nak(j),' fourth=',i

!---------------------------------------------------------------------

! (To skip over the output statements, un-comment the IF---ENDIF)


            IF (i .LT. 1) THEN
               WRITE(87,*) '#'
               WRITE(87,*) '# i=',i,' /',maxiter
               WRITE(87,*) '# state=',jstate,'  ',np(j),nh(j),nak(j)
               WRITE(87,*) '# r    ta    pdf    qdf    ph    qh'
               WRITE(87,*) '#'
               CLOSE(87)
               WRITE(npchar,'(I1)') np(j)
               WRITE(nakchar,'(I2)') nak(j)
               open(87, file=npchar//nh(j)//nakchar)
               DO lo = 1, mtp
                  WRITE(87,'(6e14.6)') r(lo), ta(lo), 
     &                           ptemp(lo),qtemp(lo),ph(lo),qh(lo)
               ENDDO
               CLOSE(87)
            ENDIF
!---------------------------------------------------------------------
*
*   Determine the first contribution to the self energy; the common
*   factor  ZEFF**4/PI*C**3 is included as the last step
*
            VALU =  OVRLAP*OVRLAP
     :             *FZALF (NPJ,KAPPA,ZEFF)
     :             /DBLE (NPJ**3)
*
*   Extract the piece that made this contribution
*
            DO 10 I = 2,MTP
               PTEMP(I) = PTEMP(I)-OVRLAP*PH(I)
               QTEMP(I) = QTEMP(I)-OVRLAP*QH(I)
   10       CONTINUE
*
*   Determine the remaining contributions to the self energy
*
            IF (KAPPA .LT. 0) THEN
               KSTART = -KAPPA
            ELSE
               KSTART = KAPPA+1
            ENDIF
            DO 12 K = KSTART,8
               IF (K .NE. NPJ) THEN
                  OVRLAP = HOVLAP (PTEMP,QTEMP,MFJ,K,KAPPA,ZEFF)
                  IF (ABS (OVRLAP) .LE. 1.0D-02) GOTO 13
                  VALU = VALU+ OVRLAP*OVRLAP
     :                        *FZALF (K,KAPPA,ZEFF)
     :                        /DBLE (K**3)
                  DO 11 I = 2,MTP
                     PTEMP(I) = PTEMP(I)-OVRLAP*PH(I)
                     QTEMP(I) = QTEMP(I)-OVRLAP*QH(I)
   11             CONTINUE
               ENDIF
   12       CONTINUE
*
*   Weight the result, multiply it by the overall factor,
*   and transfer it to the appropriate array element
*
   13       SLFINT(J) = VALU*ZEFF**4/(PI*C**3)
*
         ELSE
*
*   The self-energy for orbitals with principal quantum number
*   greater than 8 is set to zero
*
            SLFINT(J) = 0.0D 00
*
         ENDIF
*
   14 CONTINUE
*
*   Deallocate storage for the `generalised occupation numbers'
*
      CALL DALLOC (PNTUCF)
*
      RETURN
      END
