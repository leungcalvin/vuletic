************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***           **   **  *******   *****    *****    *****             ***
***           **   **  **       **   **  **   **  **   **            ***
***           **   **  **       **       **   **      **             ***
***           *******  ****      *****    *****      **              ***
***           **   **  **            **      **     **               ***
***           **   **  **       **   **     **     **                ***
***           **   **  **        *****    **      *******            ***
***                                                                  ***
***            Relativistic Hyperfine Structure Program              ***
***                         GRASP92 Version                          ***
***                         Dynamic Storage                          ***
***                                                                  ***
***   ============================================================   ***
***   Copyright (c) 1995 by P Jonsson, F A Parpia, and C F Fischer   ***
***   ============================================================   ***
***   All rights reserved.  No part of this software or its accom-   ***
***   panying documentation may be reproduced, stored in a retrie-   ***
***   val system,  or transmitted,  in any form or  by any  means,   ***
***   electronic, mechanical,  photocopying,  recording, or other-   ***
***   wise, without the prior written permission of the authors.     ***
***                                                                  ***
***                           Disclaimer                             ***
***                           ==========                             ***
***   The  authors make  no warranties,  express or implied,  that   ***
***   this software or its  accompanying documentation are free of   ***
***   error or that  they will meet your requirements for any par-   ***
***   ticular application.  Neither this software nor its accompa-   ***
***   nying documentation should be relied  upon for solving prob-   ***
***   lems if an incorrect solution could result in injury or loss   ***
***   of, or damage to, property. If you use this software or  its   ***
***   accompanying documentation,  you do so entirely  at your own   ***
***   risk;  the authors disclaim all liability for direct or con-   ***
***   sequential damage.                                             ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM HFS92
*                                                                      *
*   Entry routine for HFS92. Controls the entire computation.          *
*                                                                      *
*   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
*               [HFS92]: CHKPLT, GETHFD, HFS, SETDBG, SETSUM,          *
*                        STRSUM.                                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN, YES
      CHARACTER*24 NAME
      COMMON/DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde

      WRITE (istde,*)
      WRITE (istde,*) 'HFS92: Execution begins ...'

      WRITE (istde,*)
      WRITE (istde,*) 'Default settings?'
      YES = GETYN ()
      WRITE (istde,*)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 WRITE (istde,*) 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         WRITE (istde,*) 'Names may not start with a blank'
         GOTO 10
      ENDIF
      WRITE (istde,*)
      WRITE (istde,*) 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
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
      CALL SETSUM(NAME,NCI)
*
*   Open, check, load data from, and close, the  .csl  file
*
      CALL SETCSLA(NAME,ncore_not_used)
*
*   Get the remaining information
*
      CALL GETHFD(NAME)
*
*   Get the eigenvectors
*
      IF (NCI.EQ.0) THEN
         CALL GETMIXC(NAME)
      ELSE
         CALL GETMIXA(NAME)
      ENDIF
*
*   Append a summary of the inputs to the  .sum  file
*
      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the HFS calculation
*
      CALL HFS
*
*   Print completion message
*
      WRITE (istde,*)
      WRITE (istde,*) 'HFS92: Execution complete.'
*
      STOP
      END