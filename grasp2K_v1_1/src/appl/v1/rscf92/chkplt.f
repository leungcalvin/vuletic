************************************************************************
*                                                                      *
      SUBROUTINE CHKPLT
*                                                                      *
*   This code checks for  consistent substitutions of plants between   *
*   RSCF92 and the LIB92 subprograms.                                  *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, LODPLT.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 09 Dec 1992   *
*                                                                      *
************************************************************************
*

!$Id: chkplt.f,v 1.1 2003/09/30 05:51:39 georgio Exp $
!$Log: chkplt.f,v $
!Revision 1.1  2003/09/30 05:51:39  georgio
!
!added
!
!Revision 1.2  1997/03/04 00:15:50  xhh
!*** empty log message ***
!
      CHARACTER*(*)    RCSID
      PARAMETER        ( RCSID
     & ='$Id: chkplt.f,v 1.1 2003/09/30 05:51:39 georgio Exp $'
     & )

! RCSID etc. added
! Output lines joined
! Output quantities KEYORB, NNNP, NNNW, NNNWP released from the quoted 
!   string
! XHH 1997.01.21

      LOGICAL LPLANT
      CHARACTER*256 RECORD
      include 'parameters.def'
CGG      INTEGER KEYORB
CGG      PARAMETER (KEYORB = 121)
*
CGG      INTEGER NNNP
CGG      PARAMETER (NNNP = 590)      
CGG      INTEGER NNN1
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
      COMMON/LIB92P/LPLANT,NPLANT(4)
      COMMON/iounit/istdi,istdo,istde
*
*   Load COMMON/LIB92P/
*
      CALL LODPLT
*
*   Consistently DOUBLEPRECISION or REAL?
*
      IF (LPLANT .NEQV. .TRUE.) THEN
         IF (LPLANT) THEN
            WRITE(istde,*) 'Plant DB was set to .TRUE. in LIB92,'
     &, ' but to .FALSE. in RSCF92.'
         ELSE
            WRITE(istde,*) 'Plant DB was set to .FALSE. in LIB92,'
     &, ' but to .TRUE. in RSCF92.'
         ENDIF
         STOP
      ENDIF
*
*   Consistent numerical plants?
*
      IF (NPLANT(1) .NE. KEYORB) THEN
         CALL CONVRT (NPLANT(1),RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant KEYORB has been set to '
     &, ' '//RECORD(1:LENTH)//' in LIB92, but to ', KEYORB
     &, ' in RSCF92.'
         STOP
      ENDIF
*
      IF (NPLANT(2) .NE. NNNP) THEN
         CALL CONVRT (NPLANT(2),RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant NP has been set to '
     &, ' '//RECORD(1:LENTH)//' in LIB92, but to ',NNNP
     &, ' in RSCF92.'
         STOP
      ENDIF
*
      IF (NNN1 .NE. NNNP+10) THEN
         CALL CONVRT (NNNP+10,RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant N1 should be set'
     &, ' to '//RECORD(1:LENTH)//'.'
         STOP
      ENDIF
*
      IF (NPLANT(3) .NE. NNNW) THEN
         CALL CONVRT (NPLANT(3),RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant NW has been set to '
     &, ' '//RECORD(1:LENTH)//' in LIB92, but to ',NNNW
     &, ' in RSCF92.'
         STOP
      ENDIF
*
      IF (NPLANT(4) .NE. NNNWP) THEN
         CALL CONVRT (NPLANT(4),RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant NWP has been set to '
     &, ' '//RECORD(1:LENTH)//' in LIB92, but to ',NNNWP
     &, ' in RSCF92.'
         STOP
      ENDIF
*
      IF (MOD (NNNW,4) .EQ. 0) THEN
         NWP = NNNW/4
      ELSE
         NWP = NNNW/4+1
      ENDIF
*
      IF (NNNWP .NE. NWP) THEN
         CALL CONVRT (NWP,RECORD,LENTH)
         WRITE(istde,*) 'CHKPLT: Plant NWP should be set'
     &, ' to '//RECORD(1:LENTH)//'.'
         STOP
      ENDIF
*
      RETURN
      END
