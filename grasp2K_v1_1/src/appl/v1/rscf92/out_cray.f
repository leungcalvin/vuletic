************************************************************************
*                                                                      *
      SUBROUTINE OUT (J,JP,P,Q)
*                                                                      *
*   This subroutine carries out the step-by-step outward integration   *
*   of a pair of inhomogeneous Dirac radial equations.                 *
*                                                                      *
*   arguments:                                                         *
*                                                                      *
*      J:   (Input) Orbital index of function to be computed           *
*      JP:  (Input) The join point; the outward integration stops      *
*           at this tabulation index                                   *
*      P,Q: (Input and output) on input, elements 1 to 3 of both       *
*           arrays must be tabulated; on output, the arrays are        *
*           tabulated up to point JP                                   *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION P(NNNP),Q(NNNP)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT3/TF(NNNP),TG(NNNP),XU(NNNP),XV(NNNP)
     :      /ORB4/NP(NNNW),NAK(NNNW)

! CRAY coding - Faster on CRAY but slower on IBM RS/6000
! XHH 1997.03.07

      DIMENSION aux1(NNNP), aux2(NNNP), aux3(NNNP), 
     &          aux4(NNNP), aux5(NNNP), aux6(NNNP)

      w = 0.5d0 * H * DBLE( nak(j) )

      DO i = 4, jp

         tmp1 = 1.d0/((1.d0-w*rpor(i))*(1.d0+w*rpor(i)) - tf(i)*tg(i))

         aux1(i) = tmp1 * 
     &       ( (1.d0-w*rpor(i)) * (1.d0-w*rpor(i-1)) + tf(i) * tg(i-1) )

         aux2(i) = tmp1 *
     &       ( (1.d0-w*rpor(i)) * tf(i-1) + (1.d0+w*rpor(i-1)) * tf(i) )

         aux3(i) = tmp1 *
     &       ( (1.d0-w*rpor(i)) * xu(i-1) - tf(i) * xv(i-1) )


         aux4(i) = tmp1 *
     &       ( (1.d0+w*rpor(i)) * (1.d0+w*rpor(i-1)) + tg(i) * tf(i-1) )

         aux5(i) = tmp1 *
     &       ( (1.d0+w*rpor(i)) * tg(i-1) + (1.d0-w*rpor(i-1)) * tg(i) )

         aux6(i) = tmp1 *
     &       ( (1.d0+w*rpor(i)) * xv(i-1) - tg(i) * xu(i-1) )


      ENDDO

      DO i = 4, jp
         p(i) = p(i-1) * aux1(i) - q(i-1) * aux2(i) + aux3(i)
         q(i) = q(i-1) * aux4(i) - p(i-1) * aux5(i) + aux6(i)
      ENDDO
 
      RETURN
      END

