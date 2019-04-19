
************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***        *****   ******  **   **  **    **   *****   ******        ***
***       **   **  **      ***  **  ***  ***  **   **  **   **       ***
***       **       **      ***  **  ** ** **  **       **   **       ***
***       **  ***  ****    ** ****  ** ** **  **       ******        ***
***       **   **  **      **  ***  **    **  **       **            ***
***       **   **  **      **   **  **    **  **   **  **            ***
***        *****   ******  **   **  **    **   *****   **            ***
***                                                                  ***
***   Program for generating the energy expression for H(DC). This   ***
***   program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,  ***
***   and C. F. Fischer, 1990).                                      ***
***                                                                  ***
***                         GRASP92 Version                          ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM GENMCP
*                                                                      *
*   Entry routine for GENMCP. Controls the entire computation.         *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC.                                        *
*               [GENMCP]: CHKPLT, MATRIX, MCP, SETCSL, SETDBG,         *
*                         SETMCP, SETSUM, STRSUM.                      *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde

      LOGICAL DEBUG,RESTRT,GETYN,YES
      COMMON/DEFAULT/NDEF
*
      PRINT *, 'GENMCP: Execution begins ...'
      PRINT *

      WRITE(istde,*) 'Default settings?'
      YES = GETYN ()
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF
*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT ('GENMCP')
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG (DEBUG)
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Open the  .sum  file
*
      IF (NDEF.NE.0) CALL SETSUM
*
*   Open, check, load data from, and close the  .csl  file
*
      CALL SETCSL(NCORE_NOT_USED)
*
*   Set up the  .mcp  files; determine if this is a restart
*
      CALL SETMCP (RESTRT)
*
*   Append a summary of the inputs to the  .sum  file
*
      IF (NDEF.NE.0) CALL STRSUM
*
*   Set up the table of logarithms of factorials for use by
*   angular modules
*
      CALL FACTT
*
*   Proceed with the generation of MCP coefficients
*
      CALL MCP (RESTRT)
*
*   All done; close files that are still open
*
      CLOSE (24)
      IF (DEBUG) CLOSE (99)
*
*   Print completion message
*
      PRINT *, 'GENMCP: Execution complete.'
*
      STOP
      END
