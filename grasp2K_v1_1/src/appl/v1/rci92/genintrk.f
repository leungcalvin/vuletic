
************************************************************************
*                                                                      *
        SUBROUTINE GENINTRK
*                                                                      *
*       Generate the list of Rk Integrals that could arise from a set  *
*       of orbitals.                                                   *
*                                                                      *
*     Written by Per Jonsson                                           *
*                                                                      *
************************************************************************

        IMPLICIT DOUBLEPRECISION (A-H, O-Z)

        include 'parameters.def'
CGG        PARAMETER (NNNW = 120) 
        PARAMETER (KMAX = 20)
C        PARAMETER (KEY = 121)

        LOGICAL GEN,TRIANGRK
        POINTER(PNTRIQ,RIQDUMMY(1))

        POINTER (PCTEVLRK,VALTEIRK(1))
        POINTER (PCTEILRK, INDTEIRK(1))

        COMMON/ORB2/NCF,NW,PNTRIQ
     :        /ORB4/NP(NNNW),NAK(NNNW)
     :        /ORB5/NKL(NNNW),NKJ(NNNW)  
     :        /CTEILSRK/PCTEILRK,PCTEVLRK
     :        /KKSTART/KSTART(0:KMAX)   

        KEY = NW + 1
        KSTART(0) = 1
*
*   Sweep through to find dimensions
*
        WRITE(*,*) 'NW',NW
	GEN = .FALSE.
*
*   Find 2*JMAX
*
        J2MAX = NKJ(1)
        DO I = 2,NW
          IF (NKJ(I).GT.J2MAX) J2MAX = NKJ(I)
        ENDDO

        IF (J2MAX.GT.KMAX) THEN
          WRITE(*,*) 'The subroutine genrkint must be changed'
          STOP
        ENDIF
*
*   Make the RK integrals: IA <= IB, IA <= IC, IA <= ID, IB <= ID
*
  999   N = 0
        DO K = 0,J2MAX
          DO IA = 1,NW
            DO IB = IA,NW
              DO IC = IA,NW
                IF (TRIANGRK(NKL(IA),K,NKL(IC))) THEN
                  DO ID = IB,NW
                    IF (TRIANGRK(NKL(IB),K,NKL(ID))) THEN
                      N = N + 1
                      IF (GEN) THEN
                        INDTEIRK(N) = ((IA*KEY+IB)*KEY+IC)*KEY+ID
                        VALTEIRK(N) = SLATER(IA,IB,IC,ID,K)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          KSTART(K+1) = N + 1
        ENDDO
*
*     Allocate memory for integral book keeping
*
        IF (.NOT. GEN) THEN
          CALL ALLOC (PCTEILRK,N,4)
          CALL ALLOC (PCTEVLRK,N,8)
          WRITE(*,*) 'Allocating space for ',N,' Rk integrals'
          GEN = .TRUE.
	  GOTO 999
	END IF
*
        RETURN
        END
