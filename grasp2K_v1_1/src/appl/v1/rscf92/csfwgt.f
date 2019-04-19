************************************************************************
      SUBROUTINE CSFWGT (LSTDIO)
*                                                                      *
*   Print  the  weights of the largest five CSF contributors to each   *
*   ASF.                                                               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, ISPAR, ITJPO.          *
*                                                                      *
*                                          Last updated: 21 Dec 1992   *
*                                          Last updated: 24 Feb 1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

!A bug on nelt fixed (see 97.04.04)

!$Id: csfwgt.f,v 1.1 2003/09/30 05:51:39 georgio Exp $
!$Log: csfwgt.f,v $
!Revision 1.1  2003/09/30 05:51:39  georgio
!
!added
!
!Revision 1.3  1997/06/02 21:58:33  xhh
!*** empty log message ***
!
!Revision 1.2  1997/03/04 00:15:50  xhh
!new algorithm for 5 largest components searching
!
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LSTDIO
      CHARACTER*256 RECORD
      CHARACTER*8 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
*
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))

      COMMON/DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVEC/PNEVEC
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /SYMA/PIATJP,PIASPA

      DIMENSION coeff(5),ICONF(5)
      COMMON/iounit/istdi,istdo,istde

      IF(lstdio) THEN
         WRITE(istdo,300)
      ELSE
         WRITE(24,300)
      ENDIF

      NELT = MIN(5,NCF)
!      print *, 'ncf=',ncf, '   ncmin=',ncmin,'  nelt=',nelt

      DO 40 IV = 1, NCMIN
!         loop over eigenvectors

         ICF = ICCMIN(IV)
         ivtjpo = iatjpo(iv)    ! j-value related
         ivspar = iaspar(iv)    ! parity related
         ii = (iv-1)*ncf

         icount = 0
         DO I = 1, ncf
            IF(itjpo(i) .EQ. ivtjpo .AND. ispar(i) .EQ. ivspar) THEN
               icount = icount + 1
               coeff(icount) = evec(ii+i)
               iconf(icount) = i
               IF(icount .EQ. nelt) GOTO 10
!               normal exit
            ENDIF
         ENDDO
!         abnormal exit (should not happen)
!         normal ! If the number of CSF's for a given parity and J is
!         smaller than nelt. This happens in a mixed calculation.
!         Thus the following write and stop statements commented _and_
!         nelt has to be re-assigned with this smaller icount.
!         nelt may now be different for dofferent (parity,J) subspace
!         in a mixed job.

          nelt = icount

!XHH 97.04.04
!         WRITE(istde,*) 'csfwgt: nelt _seems_ too big'
!         STOP
   10    CONTINUE
         igroup = i

!         sort the first nelt in decreasing order
         Do i = 1, nelt
         Do j = i+1, nelt
            if (abs(coeff(j)) .gt. abs(coeff(i))) then
               temp = coeff(i)
               coeff(i) = coeff(j)
               coeff(j) = temp
               itemp = iconf(i)
               iconf(i) = iconf(j)
               iconf(j) = itemp
            end if
         end do
         end do

         Do i = igroup+1 , ncf
            IF(itjpo(i) .NE. ivtjpo .OR. ispar(i) .NE. ivspar) GOTO 20
            w = evec(ii+i)
            if (w .ne. 0.d0 .and. abs(w) .gt. abs(coeff(nelt)) ) then
!               we have a non-zero value larger than the largest so far
               do j = 1, nelt
                  if (abs(w) .gt. abs(coeff(j))) then
                     do irem = nelt, j+1, -1
                        coeff(irem) = coeff(irem-1)
                        iconf(irem) = iconf(irem-1)
                     enddo
                     coeff(j) = w
                     iconf(j) = i
                     goto 20
                  endif
               enddo
            endif
   20       CONTINUE
         end do

        ip = (iaspar(iv) + 3) / 2

         IF (LSTDIO) THEN
            RECORD(1:7) = ' J P = '
            IBEG = 8
            IEND = IBEG+6
            RECORD(IBEG:IEND) = LABJ(IATJPO(IV))//LABP(IP)(1:2)//':'
            DO 30 I = 1,NELT
               IBEG = IEND+1
               IEND = IBEG
               RECORD(IBEG:IBEG) = ' '
               IBEG = IEND+1
               IEND = IBEG+5 !was +5
               WRITE (RECORD(IBEG:IEND),'(F6.4)') coeff(I)
               IBEG = IEND+1
               IEND = IBEG+3
               RECORD(IBEG:IEND) = ' of '
               CALL CONVRT (ICONF(I),CNUM,LENTH)
               IBEG = IEND+1
               IEND = IBEG+LENTH-1
               RECORD(IBEG:IEND) = CNUM(1:LENTH)
   30       CONTINUE
            WRITE (istdo,310) RECORD(1:14), RECORD(16:IEND)
         ELSE
         
            WRITE (24,320) ICF,LABJ(IATJPO(IV)),LABP(IP),
     :                     (coeff(I),ICONF(I),I = 1,NELT)
         ENDIF
   40 CONTINUE

      RETURN

  300 FORMAT (/'Weights of major contributors to ASF:'
     :       //'Level  J Parity      CSF contributions'/)
  310 FORMAT (1X, A14,80A)
  320 FORMAT (1X,I4,4X,2A4,5(1X,F6.4,' of ',I5))

      END

