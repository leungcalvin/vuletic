************************************************************************
*                                                                      *
      SUBROUTINE dsuvec (EOL,NCF)
!-------------------------------------------------
!
! This subroutine does the following
!
! DO IR = 1, NCF
!   dsubrss(IR) = dsubrs (EOL,IR,IR)
! ENDDO
!
!-------------------------------------------------
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      POINTER (PCCMIN,CCMINDUMMY)
      POINTER (PWEIGH,WEIGHDUMMY)
      LOGICAL EOL
*
      POINTER (PNTRWT,WT(1))
      POINTER (PNEVEC,EVEC(1))
*
      COMMON/DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVEC/PNEVEC

      POINTER (pnlocal,ddsubrs(1))
      COMMON /dsu2cof/pnlocal

      IF (EOL) THEN

         DO ir = 1, NCF
            ddsubrs(ir) = 0.D0
         ENDDO

         DO  K = 1, NCMIN
            wtk = wt(k)
            kn = (K-1)*NCF
            DO ir = 1, NCF
               ddsubrs(ir) = ddsubrs(ir) + EVEC(ir+kn)*EVEC(ir+kn)*wtk
            ENDDO
         ENDDO

      ELSE

         DO ir = 1, NCF
            ddsubrs(ir) = WT(ir)
         ENDDO

      ENDIF
*
      RETURN
      END
