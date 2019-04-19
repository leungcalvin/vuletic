************************************************************************
*                                                                      *
      SUBROUTINE GETSCD (EOL)
*                                                                      *
*   Interactively determines the data governing the SCF problem.       *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, GETRSL, GETYN, NUCPOT, RADGRD,        *
*                        SETISO, SETQIC, SETRWF.                       *
*               [RSCF92]: GETALD, GETOLD, PRTRSL.                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Dec 1992   *
*                                                                      *
************************************************************************
! Short output lines joined
! PRINT *, --> WRITE(istde,*), where istde is defined in io_units
! XHH 1997.01.22

      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)

      COMMON/iounit/istdi,istdo,istde

Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG,EOL,GETYN,LFIX,LFORDR,NOINVT,ORTHST,YES
      CHARACTER*80 RECORD
      CHARACTER*20 CNUM
      CHARACTER*2 NH
*
      DIMENSION INDEX(NNNW)
*
      POINTER (PCDAMP,CDAMP(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
! PCCMIN is not used, but ncmin is.
      POINTER (PCCMIN,ICCMIN(1))
      COMMON/COUN/THRESH
     :      /DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INVT/NOINVT(NNNW)
     :      /MCPB/DIAG,LFORDR
      COMMON/NODE/NNODEP(NNNW)
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORBA/IORDER(NNNW)
     :      /ORTHCT/ORTHST
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Open, check, load data from, and close the  .iso  file
*
      CALL SETISO('isodata')
*
*   Set default speed of light and grid parameters
*
      C = CVAC
      IF (NPARM .EQ. 0) THEN
         RNT = EXP (-65.0D 00/16.0D 00) / Z
         H = 0.5D 00**4
         N = MIN (220,NNNP)
      ELSE
!         default comes here
         RNT = 2.0D-06
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D 00
*
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'Change the default speed of light'
     &, ' or radial grid parameters?'
         YES = GETYN ()
*
         IF (YES) THEN
*
*   Revise the speed of light
*
            WRITE(istde,*) 'The physical speed of light in'
     &, ' atomic units is',CVAC,';'
     &, ' revise this value?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE(istde,*) 'Enter the revised value:'
               READ *,C
            ENDIF
*
*   Determine the parameters controlling the radial grid
*
            WRITE(istde,*) 'The default radial grid parameters for '
     &,                    'this case are:'
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
*   Load the subshell radial wavefunction estimates
*
      CALL SETRWF
*
*   Set some defaults
*
      THRESH = 0.05D 00
*
      DO 1 I = 1,NW
!XHH       The following 5 are to be modified in  getold 
!          Note initializations for METHOD, NOINVT and ODAMP have been
!          changed due to the change of the question from 
!              'Which of these are correlation functions?'
!          to   
!              'Which of these are NOT correlation functions?'
!          in subroutine getold.f
!
         IORDER(I) = I
         METHOD(I) = 3
         NOINVT(I) = .TRUE.
         ODAMP(I) = 1.0D 00
         SCNSTY(I) = 0.0D 00

         IF (NAK(I) .LT. 0) THEN
            NNODEP(I) = NP(I)+NAK(I)
         ELSE
            NNODEP(I) = NP(I)-NAK(I)-1
         ENDIF
    1 CONTINUE
*
*
      IF (DIAG) THEN
         EOL = .FALSE.
         WRITE(istde,*) '(E)AL type calculation; H(DC) will not be '
     &,                 'diagonalised;'
         CALL GETALD
      ELSEIF (LFORDR) THEN
         EOL = .TRUE.
         WRITE(istde,*) '(E)OL type calculation;'
         CALL GETOLD
      ELSE
         WRITE(istde,*) '(E)OL type calculation?'
         YES = GETYN ()
         IF (YES) THEN
            EOL = .TRUE.
            CALL GETOLD
         ELSE
            EOL = .FALSE.
            WRITE(istde,*) '(E)AL type calculation; H(DC) will not be '
     &,                    'diagonalised;'
            CALL GETALD
         ENDIF
      ENDIF

!      CALL CONVRT (NSCF,CNUM,LENTH)
      WRITE(istde,*) 'Enter the maximum number of SCF cycles:'
      READ (*,*) NSCF
!         CALL updord
*
*   Allow the user to modify other defaults
*
      IF (NDEF.NE.0) THEN
         WRITE(istde,*) 'Modify other defaults?'
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF
*
      IF (YES) THEN
*
*   THRESH
*
         WRITE(istde,*) 'An oscillation in the large-component of the '
     &,                 'radial wavefunction is diregarded'
         WRITE(istde,*) 'for the purposes of node counting if its '
     &,                 'amplitude is less than 1/20 the'
         WRITE(istde,*) 'maximum amplitude.   Revise this?'
         YES = GETYN ()
         IF (YES) THEN
    3       WRITE(istde,*) 'Enter the new threshold value:'
            READ (*,*) THRESH
            IF (THRESH .LE. 0.0D 00) THEN
               WRITE(istde,*) 'GETSCD: This must exceed 0;'
               GOTO 3
            ENDIF
         ENDIF
*
*   IORDER - has been determined at the place where the varying 
*            orbitals are specified
*
!         WRITE(istde,*) 'The subshells will be improved in the order'
!         CALL PRTRSL
!         WRITE(istde,*) 'Revise this order?'
!         YES = GETYN ()
         YES = .FALSE.
         IF (YES) THEN
            WRITE(istde,*) 'Revised order:'
    4       CALL GETRSL (INDEX,NSUBS)
            NORDER = 0
            DO 6 I = 1,NSUBS
               NPLOC = NP(INDEX(I))
               NAKLOC = NAK(INDEX(I))
               DO 5 J = 1,NW
                  IF ((NP(J) .EQ. NPLOC) .AND.
     :                (NAK(J) .EQ. NAKLOC) .AND.
     :                (.NOT. LFIX(J))) THEN
                     NORDER = NORDER+1
                  ENDIF
    5          CONTINUE
    6       CONTINUE
            IF (NORDER .NE. NW-NFIX) THEN
               WRITE(istde,*) 'GETSCD: All subshells that are to be '
     &,                       'improved must appear in the list.'
               GOTO 4
            ELSE
               NORDER = 0
               DO 7 I = 1,NW
                  IF (.NOT. LFIX(I)) THEN
                     NORDER = NORDER+1
                     IORDER(I) = INDEX(NORDER)
                  ENDIF
    7          CONTINUE
            ENDIF
         ENDIF
*
*   METHOD
*
         DO 11 I = 1,4
            DO 8 J = 1,NW
               IF ((METHOD(J) .EQ. I) .AND. (.NOT. LFIX(J))) THEN
!                  CALL CONVRT (I,CNUM,LENTH)
                  WRITE(istde,*) 'Method ',I,' is used for '
     &,                          'integrating the radial differential '
     &,                          'equation for subshells'
                  GOTO 9
               ENDIF
    8       CONTINUE
            GOTO 11
    9       IEND = 0
            DO 10 J = 1,NW
               IF ((METHOD(J) .EQ. I) .AND. (.NOT. LFIX(J))) THEN
                  IBEG = IEND+1
                  IEND = IBEG
                  RECORD(IBEG:IEND) = ' '
                  CALL CONVRT(NP(J),CNUM,LENTH)
                  IBEG = IEND+1
                  IEND = IBEG+LENTH-1
                  RECORD(IBEG:IEND) = CNUM(1:LENTH)
                  IBEG = IEND+1
                  IF (NAK(J) .LT. 0) THEN
                     IEND = IBEG
                     RECORD(IBEG:IEND) = NH(J)(1:1)
                  ELSE
                     IEND = IBEG+1
                     RECORD(IBEG:IEND) = NH(J)(1:2)
                  ENDIF
               ENDIF
               IF (IEND .GT. 76) THEN
                  WRITE(istde,*) RECORD(1:IEND)
                  IEND = 0
               ENDIF
   10       CONTINUE
            IF (IEND .GT. 0) WRITE(istde,*) RECORD(1:IEND)
   11    CONTINUE
         WRITE(istde,*) 'Select a different integration method for '
     &,                 'any subshell radial wavefunction?'
         YES = GETYN ()
         IF (YES) THEN
            DO 13 I = 1,4
!               CALL CONVRT (I,CNUM,LENTH)
               WRITE(istde,*) 'Method ',I,':'
               CALL GETRSL (INDEX,NSUBS)
               DO 12 J = 1,NSUBS
                  LOC = INDEX(J)
                  IF (.NOT. LFIX(LOC)) METHOD(LOC) = I
   12          CONTINUE
   13       CONTINUE
         ENDIF
*
*   NOINVT
*
         WRITE(istde,*) 'The first oscillation of the large component'
         DO 14 I = 1,NW
            IF (NOINVT(I) .AND. (.NOT. LFIX(I))) THEN
               WRITE(istde,*) 'of the following radial wavefunctions '
     &,                       'will be required to be positive'
               GOTO 15
            ENDIF
   14    CONTINUE
         WRITE(istde,*) 'of all radial wavefunctions will be required '
     &,                 'to be positive.   Revise this?'
         YES = GETYN ()
         GOTO 17
   15    IEND = 0
         DO 16 I = 1,NW
            IF (NOINVT(I) .AND. (.NOT. LFIX(I))) THEN
               IBEG = IEND+1
               IEND = IBEG
               RECORD(IBEG:IEND) = ' '
               CALL CONVRT(NP(I),CNUM,LENTH)
               IBEG = IEND+1
               IEND = IBEG+LENTH-1
               RECORD(IBEG:IEND) = CNUM(1:LENTH)
               IBEG = IEND+1
               IF (NAK(I) .LT. 0) THEN
                  IEND = IBEG
                  RECORD(IBEG:IEND) = NH(I)(1:1)
               ELSE
                  IEND = IBEG+1
                  RECORD(IBEG:IEND) = NH(I)(1:2)
               ENDIF
            ENDIF
            IF (IEND .GT. 76) THEN
               WRITE(istde,*) RECORD(1:IEND)
               IEND = 0
            ENDIF
   16    CONTINUE
         IF (IEND .GT. 0) WRITE(istde,*) RECORD(1:IEND)
         WRITE(istde,*) 'Revise this?'
         YES = GETYN ()
   17    IF (YES) THEN
            WRITE(istde,*) 'Suppressing enforcement of positive first '
     &,                    'oscillation:'
            CALL GETRSL (INDEX,NSUBS)
            DO 18 I = 1,NSUBS
               LOC = INDEX(I)
               IF (.NOT. LFIX(LOC)) NOINVT(LOC) = .TRUE.
   18       CONTINUE
         ENDIF
*
*   ODAMP
*
         DO 19 I = 1,NW
            IF ((ODAMP(I) .NE. 0.0D 00) .AND.
     :          (.NOT. LFIX(I))) THEN
               WRITE(istde,*) 'Subshell accelerating parameters have '
     &,                       'been set.   Revise these?'
               YES = GETYN ()
               GOTO 20
            ENDIF
   19    CONTINUE
         WRITE(istde,*) 'Set accelerating parameters for subshell '
     &,                 'radial wavefunctions?'
         YES = GETYN ()
   20    IF (YES) THEN
            WRITE(istde,*) 'Different accelerating parameters for '
     &,                    'different subshell radial wavefunction?'
            YES = GETYN ()
            IF (YES) THEN
   21          WRITE(istde,*) 'Enter an accelerating parameter'
               WRITE(istde,*) ' (0< ODAMP < 1 allows ODAMP to be '
     &,                       'reduced as convergence is approached;'
               WRITE(istde,*) ' -1 < ODAMP < 0 implies |ODAMP| is '
     &,                       'held constant):'
               READ (*,*) ODAMPU
               IF ((ABS (ODAMPU) .EQ. 0.0D 00) .OR.
     :             (ABS (ODAMPU) .GE. 1.0D 00)) THEN
                  WRITE(istde,*) 'GETSCD: Value out of range ...'
                  GOTO 21
               ELSE
                  CALL GETRSL (INDEX,NSUBS)
                  DO 22 I = 1,NSUBS
                     LOC = INDEX(I)
                     IF (.NOT. LFIX(LOC)) ODAMP(LOC) = ODAMPU
   22             CONTINUE
               ENDIF
            ELSE
   23          WRITE(istde,*) 'Enter the accelerating parameter'
               WRITE(istde,*) ' (0< ODAMP < 1 allows ODAMP to be '
     &,                       'reduced as convergence is approached;'
               WRITE(istde,*) ' -1 < ODAMP < 0 implies |ODAMP| is '
     &,                       'held constant):'
               READ (*,*) ODAMPU
               IF ((ABS (ODAMPU) .EQ. 0.0D 00) .OR.
     :             (ABS (ODAMPU) .GE. 1.0D 00)) THEN
                  WRITE(istde,*) 'GETSCD: Value out of range ...'
                  GOTO 23
               ELSE
                  DO 24 I = 1,NW
                     IF (.NOT. LFIX(I)) ODAMP(I) = ODAMPU
   24             CONTINUE
               ENDIF
            ENDIF
         ENDIF
*
*   CDAMP
*
         WRITE(istde,*) 'Set accelerating parameters for the '
     &,                 'eigenvectors?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Different accelerating parameters '
     &,                    'for each eigenvector?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE(istde,*) 'Enter an accelerating parameter for'
               CALL CONVRT (NCMIN,RECORD,LENTH)
               WRITE(istde,*) ' each of the '//RECORD(1:LENTH)//
     &                        ' levels :'
               READ (*,*) (CDAMP(I),I = 1,NCMIN)
            ELSE
               WRITE(istde,*) 'Enter the accelerating parameter:'
               READ (*,*) CDAMPU
               DO 25 I = 1,NCMIN
                  CDAMP(I) = CDAMPU
   25          CONTINUE
            ENDIF
         ENDIF
*
*   NSIC
*
!         CALL CONVRT (NSIC,CNUM,LENTH)
         WRITE(istde,*) 'Following the improvement of each of the '
     &,                 'subshell radial wavefunctions in turn, '
         WRITE(istde,*) 'the ',NSIC,' least self-consistent'
     &,                 ' functions will be improved at the'
         WRITE(istde,*) 'end of the first SCF cycle. Revise this '
     &,                 'setting?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter the number of additional '
     &,                    'improvements:'
            READ (*,*) NSIC
         ENDIF
*
*   NSOLV
*
!         CALL CONVRT (NSOLV,CNUM,LENTH)
         WRITE(istde,*) 'The maximum number of cycles in attempting '
     &,                 'to solve each radial equation is '
         WRITE(istde,*) NSOLV,' times the principal quantum'
     &,                ' number of the radial'
         WRITE(istde,*) 'wave-function to be estimated.   '
     &,                 'Revise this setting?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter the factor that multiplies the '
     &,                    'principal quantum number:'
            READ (*,*) NSOLV
         ENDIF
*
*   Orthogonalisation
*
         IF (ORTHST) THEN
            WRITE(istde,*) 'Subshell radial wavefunctions will be '
     &,                    'Schmidt orthogonalised immediately'
            WRITE(istde,*) 'following their estimation to all '
     &,                    'functions with poorer self-consistency.'
            WRITE(istde,*) ' Revise this?'
            YES = GETYN ()
            IF (YES) ORTHST = .FALSE.
         ELSE
            WRITE(istde,*) 'Subshell radial wavefunctions will be '
     &,                    'Schmidt orthogonalised at the end of'
            WRITE(istde,*) 'each SCF cycle.   Revise this?'
            YES = GETYN ()
            IF (YES) ORTHST = .TRUE.
         ENDIF
*
      ENDIF
*
      RETURN
      END
