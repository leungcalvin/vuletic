************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***          *****   *******  **   **  ****   *****    *****         ***
***         **   **  **       ***  **   **   **   **  **   **        ***
***         **       **       ** * **   **   **       **   **        ***
***         **  ***  ****     ** * **   **    *****   **   **        ***
***         **   **  **       ** * **   **        **  **   **        ***
***         **   **  **       **  ***   **   **   **  **   **        ***
***          *****   *******  **   **  ****   *****    *****         ***
***                                                                  ***
***   Package to generate the nuclear charge, geometry, mass, spin,  ***
***   and electromagnetic moment data file for the GRASP92 codes.    ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM GENISO
*                                                                      *
*   Generates the isotope data file for the GRASP92 suite of codes.    *
*                                                                      *
*   Call(s) to: GETCPR, GETYN, LENGTH, OPENFL.                         *
*                                                                      *
*   Written by Farid A. Parpia.           Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      LOGICAL GETYN,YES
*
!XHH To be more implementation-independent, add 3 lines to each main 
!XHH program
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde

cb alpha constant from lib/lib92/setcon.f
cb AUMAMU from lib/lib92/setcon.f
cb
cb    DATA EMEAMU /5.48579903D-04/
cb   :     ALFAI  /137.0359895D 00/
      COMMON/DEF2/C
     :      /DEF11/FMTOAU,AUMAMU
cb
      ALFAI = 1.0D 00/C
cb EMEAMU: Electron mass in amu
      AUMAMU = EMEAMU
cb end alpha constant
*
*   File  grasp92.iso  is FORMATTED
*
      DEFNAM = 'isodata'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      FILNAM = DEFNAM
*
      CALL OPENFL (22,FILNAM,FORM,STATUS,IERR)
*
      IF (IERR .NE. 0) THEN
         WRITE(istde,*) 'Error when opening isodata'
         STOP
      ENDIF
*
      WRITE(istde,*) 'Enter the atomic number:'

      READ *, Z
      WRITE (22,300) 'Atomic number:'
      WRITE (22,*) Z
*
      WRITE(istde,*) 'Enter the mass number (0 if the'
     &, ' nucleus is to be modelled as a point source:'

      READ *, A

      WRITE (22,300) 'Mass number (integer) :'
      WRITE (22,*) A

      IF (A .EQ. 0.0D 00) THEN

         CPARM = 0.0D 00
         APARM = 0.0D 00
      ELSE
         RRMS =  0.836D 00*A**(1.0D 00/3.0D 00)
     :          +0.570D 00
         WRITE(istde,*) 'The default root mean squared'
     &, ' radius is ',RRMS,' fm;'
         TPARM = 2.30D 00
         WRITE(istde,*) ' the default nuclear skin thickness'
     &, ' is ',TPARM,' fm;'
         WRITE(istde,*) 'Revise these values?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter the root mean squared'
     &, ' radius of the nucleus (in fm):'
            READ *, RRMS
            WRITE(istde,*) 'Enter the skin thickness of'
     &, ' the nucleus (in fm):'
            READ *, TPARM
         ENDIF
         APARM = TPARM/(4.0D 00*LOG (3.0D 00))
         CALL GETCPR (RRMS,APARM,CPARM)
      ENDIF
      WRITE (22,300) 'Fermi distribution parameter a:'
      WRITE (22,*) APARM
      WRITE (22,300) 'Fermi distribution parameter c:'
      WRITE (22,*) CPARM
*
      WRITE(istde,*) 'Enter the mass of the neutral'
     &, ' atom (in amu) (0 if the nucleus is to be static):'
      READ *, AMAMU
      IF (AMAMU .NE. 0.0D 00) THEN
!         WRITE(istde,*) 'Enter your best estimate of the ground'
!     &, ' state energy of, the neutral atom (in Hartrees):'
!         READ *, EBIND
          EBIND = 0.D0

!XHH better use NINT, not INT. 
!         NENEU = INT (Z)
         NENEU = NINT (Z)

!         WRITE(istde,*) 'The number of electrons in the'
!     &, ' neutral atom is deduced to be',NENEU,';'
!         WRITE(istde,*) 'Revise this?'
!         YES = GETYN ()
!         IF (YES) THEN
!            WRITE(istde,*) 'Enter the number of electrons'
!     &, ' in the neutral atom:'
!            READ *, NENEU
!         ENDIF
         IF (EBIND .GT. 0.0D 00) EBIND = -EBIND
         EMNAMU = AMAMU-EMEAMU*DBLE (NENEU)
     :                 -EMEAMU*EBIND/ALFAI**2
      ELSE
         EMNAMU = 0.0D 00
      ENDIF
*
      WRITE (22,300) 'Mass of nucleus (in amu):'
      WRITE (22,*) EMNAMU
*
      WRITE(istde,*) 'Enter the nuclear spin quantum'
     &, ' number (I) (in units of h / 2 pi):'
      READ *, SQN
      WRITE (22,300) 'Nuclear spin (I) (in units of h / 2 pi):'
      WRITE (22,*) SQN
*
      WRITE(istde,*) 'Enter the nuclear dipole moment'
     &, ' (in nuclear magnetons):'
      READ *, DMOMNM
      WRITE (22,300) 'Nuclear dipole moment (in nuclear magnetons):'
      WRITE (22,*) DMOMNM
*
      WRITE(istde,*) 'Enter the nuclear quadrupole'
     &, ' moment (in barns):'
      READ *, QMOMB
      WRITE (22,300) 'Nuclear quadrupole moment (in barns):'
      WRITE (22,*) QMOMB
*
      CLOSE (22)
*
      STOP
*
  300 FORMAT (A)
*
      END
