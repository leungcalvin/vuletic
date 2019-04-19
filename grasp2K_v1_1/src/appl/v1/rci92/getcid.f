************************************************************************
*                                                                      *
      SUBROUTINE GETCID(NAME)
*                                                                      *
*   Interactively determines the data governing the CI problem.        *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
*               [RCI92]: SETISO, SETRWF.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*24 NAME,rwffile
      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
      COMMON/iounit/istdi,istdo,istde
      print *, 'Entering getcid...'
      rwffile=trim(name)//'.w'
*
*   Open, check, load data from, and close the  .iso  file
*
      CALL SETISO('isodata')
      print *, 'after SETISO   ...'
*
*   Determine the physical effects specifications
*
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'Revise the physical speed of light (',CVAC,
     &                  ' in a.u.) ?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter the revised value:'
            READ *,C
         ELSE
            C = CVAC
         ENDIF
      ELSE
         C = CVAC
      ENDIF
*
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'Treat contributions of some CSFs'
     &,                 ' as first-order perturbations?'
         LFORDR = GETYN ()
      ELSE
         LFORDR = .FALSE.
      ENDIF
      IF (LFORDR) THEN
         WRITE(istde,*) 'The contribution of CSFs'
         WRITE(istde,*) ' 1 -- ICCUT will be treated'
         WRITE(istde,*) ' variationally; the remainder'
         WRITE(istde,*) ' perturbatively; enter ICCUT:'
         READ *, ICCUT
      ELSE
         ICCUT = 0
      ENDIF
*
      IF (IPRERUN.EQ.0) THEN
        WRITE(istde,*) ' Do you want to prerun with limited'
     &,                ' interaction?'
        YES = GETYN ()
        IF (YES) THEN
          IPRERUN = 1
          LTRANS = .FALSE.
          LVP = .FALSE.
          LNMS = .FALSE.
          LSMS = .FALSE.
          LSE = .FALSE.
          WRITE(istde,*)  ' Give CSL cut'
          READ *, NCSFPRE
          WRITE(istde,*)  ' Give coefficient cut for H_0'
          READ *, COEFFCUT1
          WRITE(istde,*) ' Give coefficient cut for the transvers'
     &,                  ' interaction'
          READ *, COEFFCUT2
          GOTO 99
        ENDIF
      ENDIF

      WRITE(istde,*) 'Include contribution of H (Transverse)?'
      LTRANS = GETYN ()
      IF (LTRANS) THEN
         WRITE(istde,*) 'Modify all transverse photon frequencies?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter the scale factor:'
            READ *, WFACT
         ELSE
            WFACT = 1.0D 00
         ENDIF
      ELSE
         WFACT = 0.0D 00
      ENDIF
*
      WRITE(istde,*) 'Include contribution of H (Vacuum Polarisation)?'
      LVP = GETYN ()
*
      WRITE(istde,*) 'Include contribution of H (Normal Mass Shift)?'
      LNMS = GETYN ()
*
      WRITE(istde,*) 'Include contribution of H (Specific Mass Shift)?'
      LSMS = GETYN ()
*
      WRITE(istde,*) 'Estimate contributions from self-energy?'
      LSE = GETYN ()
*
*   Determine the parameters controlling the radial grid
*
  99  IF (NPARM .EQ. 0) THEN
         RNT = EXP (-65.0D 00/16.0D 00) / Z
         H = 0.5D 00**4
         N = MIN (220,NNNP)
      ELSE
         RNT = 2.0D-06
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D 00
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'The default radial grid parameters'
     &,                 ' for this case are:'
         WRITE(istde,*) ' RNT = ',RNT,';'
         WRITE(istde,*) ' H = ',H,';'
         WRITE(istde,*) ' HP = ',HP,';'
         WRITE(istde,*) ' N = ',N,';'
         WRITE(istde,*) ' revise these values?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter RNT:'
            READ *, RNT
            WRITE(istde,*) 'Enter H:'
            READ *, H
            WRITE(istde,*) 'Enter HP:'
            READ *, HP
            WRITE(istde,*) 'Enter N:'
            READ *, N
         ENDIF
      ENDIF
*
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Generate the radial grid and all associated arrays
*
      CALL RADGRD
*
*   Generate $- r \times V_ (r)$
*
      CALL NUCPOT
*
*   Load the radial wavefunctions
*
      CALL SETRWFA(rwffile)
*
*   Write the basic parameters of the model electron cloud to the
*   .res  file; this is the second record on the file --- the
*   first is the header (see SUBROUTINE SETRES)
*
      WRITE (26) NELEC,NCF,NW
*
*   Write the nuclear parameters and $- r \times V_ (r)$
*   to the  .res  file; these are the third, fourth, and fifth
*   records on the file
*
      WRITE (26) Z,EMN
      WRITE (26) NPARM,(PARM(I),I = 1,NPARM)
      WRITE (26) N,(ZZ(I),I = 1,N),NNUC
*
*   Write the physical effects specification to the  .res  file;
*   this is the sixth record on the file
*
      WRITE (26) C,LFORDR,ICCUT,LTRANS,WFACT,LVP,LNMS,LSMS
*
*   Write the grid data to the  .res  file; this is the seventh
*   record on the file
*
      NP10 = N+10
      WRITE (26) RNT,H,HP,(R(I),I = 1,NP10),(RP(I),I = 1,NP10),
     :           (RPOR(I),I = 1,NP10)
*
*   Write out the interpolated radial wavefunctions; there are
*   2*NW such records; thus the total number of records written
*   at this point is 7+2*NW
*
      DO 1 J = 1,NW
         WRITE (26) E(J),GAMA(J),PZ(J),MF(J)
         WRITE (26) (PF(I,J),I = 1,MF(J)),(QF(I,J),I = 1,MF(J))
    1 CONTINUE
*
      RETURN
      END
