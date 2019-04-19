************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***             ******    *****   ****   *****    *****              ***
***             **   **  **   **   **   **   **  **   **             ***
***             **   **  **        **   **   **       **             ***
***             ******   **        **    *****       **              ***
***             **  **   **        **      **       **               ***
***             **   **  **   **   **     **      **                 ***
***             **   **   *****   ****   **      *******             ***
***                                                                  ***
***          Relativistic Configuration-Interaction Program          ***
***                                                                  ***
***   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
***   Grant, and C. F. Fischer, 1990).                               ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM RCI92
*                                                                      *
*   Entry routine for RCI92. Controls the entire computation.          *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
*                        SETRES, SETSUM, STRSUM.                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*24 NAME
      LOGICAL GETYN,YES
      COMMON/DEFAULT/NDEF
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2

!XHH Added 3 lines:
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde

      IPRERUN = 0
      PRINT *, 'RCI92: Execution begins ...'
      PRINT *
      WRITE(istde,*) 'Default settings?'
      YES = GETYN ()
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 WRITE(istde,*) 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         WRITE(istde,*) 'Names may not start with a blank'
         GOTO 10
      ENDIF
   99 WRITE(istde,*) 'RCI92: Execution begins ...'
*
*   Check compatibility of plant substitutions
*
      PRINT *, 'Calling CHKPLLT...'
      CALL CHKPLT
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      PRINT *, 'Calling SETDBG...'
      CALL SETDBG
*
*   Perform machine- and installation-dependent setup
*
      PRINT *, 'Calling SETMC...'
      CALL SETMC
*
*   Set up the physical constants
*
      PRINT *, 'Calling SETCON...'
      CALL SETCON
*
*   Open the  .sum  file
* XHH which is named as "NAME.csum" and has a status of "NEW"
*
      PRINT *, 'Calling SETSUM...'
      CALL SETSUM(NAME)
*
*   Open, check, load data from, and close the  .csl  file
*XHH which has the name "NAME.c"
*
      PRINT *, 'Calling SETCSLA...'
      CALL SETCSLA(NAME, ncore_not_used)
*
*   Set up the  .res  file; determine if this is a restart
*
      PRINT *, 'Calling SETRES...'
      CALL SETRES(NAME)
*
*   Open the  .mix  file; determine the eigenpairs required
*
      PRINT *, 'Calling SETMIX...'
      CALL SETMIX(NAME)
*
*   Append a summary of the inputs to the  .sum  file
*
      PRINT *, 'Calling STRSUM...'
      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      PRINT *, 'Calling FACTT...'
      CALL FACTT
*
*   Calculate all the needed Rk integrals
*
      PRINT *, 'Calling GENINTRK...'
      CALL GENINTRK
*
*   Proceed with the CI calculation
*
      PRINT *, 'Calling MATRIX...'
      CALL MATRIX
*
      IF (IPRERUN.EQ.1) THEN
         IPRERUN = 2
         GOTO 99
      ENDIF
*
*   Print completion message
*
      PRINT *, 'RCI92: Execution complete.'
*
      STOP
      END
