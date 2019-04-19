************************************************************************
*                                                                      *
      FUNCTION DSUBRS (EOL,I,J)
*                                                                      *
*   The coefficients d   for  I = r, J = s  are calculated here.       *
*                     rs                                               *
*                                                                      *
*                         NCMIN                                        *
*                          Sum     w     c          c                  *
*                         t = 1     t     r Gamma    s Gamma           *
*                                                t          t          *
*    (E)OL:        d   = -------------------------------------         *
*                   rs             NCMIN                               *
*                                   Sum     w                          *
*                                  t = 1     t                         *
*                                                                      *
*                                                                      *
*                                           /  NCF                     *
*    (E)AL:           d   =  delta     w   /   Sum   w                 *
*                      rs         rs    r /   t = 1   t                *
*                                                                      *
*                                                                      *
*   Written by Farid A Parpia               Last update: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
Cww      INTEGER PCCMIN,PNTRIQ,PWEIGH
      POINTER (PCCMIN,CCMINDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PWEIGH,WEIGHDUMMY)
      LOGICAL EOL
*
      POINTER (PNTRWT,WT(1))
      POINTER (PNEVEC,EVEC(1))
*
      COMMON/DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
*
      IF (EOL) THEN
*
         DSUBRS = 0.0D 00
         DO 1 K = 1,NCMIN
            DSUBRS = DSUBRS+EVEC(I+(K-1)*NCF)
     :                     *EVEC(J+(K-1)*NCF)*WT(K)
    1    CONTINUE
*
      ELSE
*
         IF (I .EQ. J) THEN
            DSUBRS = WT(I)
         ELSE
            DSUBRS = 0.0D 00
         ENDIF
*
      ENDIF
*
      RETURN
      END
