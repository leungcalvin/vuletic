************************************************************************
*                                                                      *
      LOGICAL FUNCTION TRIANGRK (LA,K,LB)
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
*   Perform the triangularity check
*
      IF (MOD(K+LA+LB,2) .NE. 0) THEN
        TRIANGRK = .FALSE.
      ELSE
        IF (ABS(LA-LB) .GT. K ) THEN
          TRIANGRK = .FALSE.
        ELSE IF (LA+LB .LT. K) THEN
          TRIANGRK = .FALSE.
        ELSE
          TRIANGRK = .TRUE.
        ENDIF
      ENDIF
 
      RETURN
      END
