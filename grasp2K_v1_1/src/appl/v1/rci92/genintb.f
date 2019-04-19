************************************************************************
*                                                                      *
      SUBROUTINE GENINTB
*                                                                      *
*     Generate the list of Breit integrals that could arise from       *
*     a set of orbitals.                                               *
*                                                                      *
*     Written by Per Jonsson                                           *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120) 

      LOGICAL GEN,TRIANGRK
Cww        INTEGER PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(1))

      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)  
     :      /STOR/KEEP(2,2)
     :      /BALLOC/N(2)
*                                                                       
*   Initialization for  BESSEL; this is done only once                  
*                                                                       
      KEEP(1,1) = 0                                               
      KEEP(1,2) = 0                                               
      KEEP(2,1) = 0                                               
      KEEP(2,2) = 0                                               

*
*   Sweep through to find dimensions
*
*   Find 2*JMAX
*
      J2MAX = NKJ(1)
      DO I = 2,NW
        IF (NKJ(I).GT.J2MAX) J2MAX = NKJ(I)
      ENDDO

*
*   Make the type 1 and 2 Breit integrals: 
*
      N(1) = 0
      N(2) = 0 
      DO IT = 1,2
        DO K = 0,J2MAX
          DO IA = 1,NW
            DO IB = 1,NW
              DO IC = 1,NW
                IF (TRIANGRK(NKL(IA),K,NKL(IB))) THEN
                  DO ID = 1,NW
                    IF (TRIANGRK(NKL(IC),K,NKL(ID))) THEN
                      N(IT) = N(IT) + 1
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*
      N(1) = N(1)/2
      N(2) = N(2)/2
      WRITE(*,*) 'Allocating space for ',N(1),' Breit integ. of type 1'
      WRITE(*,*) 'Allocating space for ',N(2),' Breit integ. of type 2'
*
      RETURN
      END
