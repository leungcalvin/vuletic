************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***       ******    *****    *****   *******   *****    *****        ***
***       **   **  **   **  **   **  **       **   **  **   **       ***
***       **   **  **       **       **       **   **       **       ***
***       ******    *****   **       ****      *****       **        ***
***       **  **        **  **       **          **       **         ***
***       **   **  **   **  **   **  **         **      **           ***
***       **   **   *****    *****   **        **      *******       ***
***                                                                  ***
***            Relativistic Self-Consistent-Field Program            ***
***                                                                  ***
***   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
***   Grant, and C. F. Fischer, 1990)                                ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM RSCF92
*                                                                      *
*   Entry routine for RSCF92. Controls the entire computation.         *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RSCF92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,      *
*                        SETRES, SETSUM, STRSUM.                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde

      LOGICAL EOL,GETYN,YES
      COMMON/DEFAULT/NDEF
      COMMON/CORE/NCORE
*
      WRITE (istde,*) 'RSCF92: Execution begins ...'
      WRITE (istde,*)
      WRITE (istde,*) 'Default settings?'
      YES = GETYN ()
      WRITE (istde,*)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF
*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Set up the physical constants
*
      CALL SETCON
*
*   Open the  .sum  file
*
      CALL SETSUM
*
*   Open, check, load data from, and close the  .csl  file
*
      CALL SETCSL(NCORE)
*
*   Open and check the  .mcp  files
*
      CALL SETMCP
*
*   Gather all remaining information and perform some
*   setup
*
      CALL GETSCD (EOL)
*
*   Append a summary of the inputs to the  .sum  file
*
      CALL STRSUM
*
*   Set up the files for the improved radial wavefunctions
*   and the mixing coefficients
*
! This will be done in orbout.f or scf.f
! To close the for the most of the run-time.
! XHH 1997.02.14
!      CALL SETRWG

      IF (EOL) CALL SETMIX
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the SCF calculation close all files except
*   the  .sum  file
*
      CALL SCF (EOL)
*
*   Complete the summary
*
      CALL ENDSUM
*
*   Print completion message
*
      WRITE(istde,*) 'RSCF92: Execution complete.'
*
      END
