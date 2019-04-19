************************************************************************
*                                                                      *
      SUBROUTINE GETOLD
*                                                                      *
*   Interactively determines the data governing OL problem.            *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, GETRSL, GETYN, RALLOC.         *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PNTRWT,RWTDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL GETYN,LFIX,NOINVT,ORTHST,YES
      CHARACTER*256 RECORD
      CHARACTER*7 FORM
      CHARACTER*3 CNUM
      CHARACTER*1 RECI
*
      DIMENSION INDEX(NNNW)
*
      POINTER (PCDAMP,CDAMP(1))
      POINTER (PWEIGH,WEIGHT(1))
      POINTER (PCCMIN,ICCMIN(1))
*
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /DEFAULT/NDEF
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /INVT/NOINVT(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORBA/IORDER(NNNW)
     :      /ORTHCT/ORTHST
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
      COMMON/iounit/istdi,istdo,istde
      LOGICAL lcorre
      COMMON/corre/lcorre(NNNW)
      SAVE /corre/
*
*   Make the initial allocation for the level list vector and
*   the auxiliary vector required by SUBROUTINE NEWCO
*
      NCD = 1
      CALL ALLOC (PCCMIN,NCD,4)
*
*   Entry message
*
    1 WRITE(istde,*) 'Enter the serial numbers of the ASF(s)'
*
*   Read and parse the list of levels; this takes a lot of code
*   because the number of levels isn't known before the list is
*   typed in
*
      READ (*,'(A)') RECORD
*
*   Initialise NCMIN
*
      NCMIN = 0
*
*   Parse RECORD from left to right
*   
      IFIRST = 0
      ISTART = 1
      I = 1
*
*   .. skip the blanks and commas(this implementation allows input to 
*      start with blanks
    2 RECI = RECORD(I:I)
      IF ((RECI .NE. ' ') .AND. (RECI .NE. ',')) THEN
	istart = i
      ELSE
	i = i+1
	if (i .le. 256) then
	  go to 2
	else
	  go to 4
	end if
      END IF
*   .. search for end of string (blank, comma, or dash)
    3 RECI = RECORD(I:I)
      IF ((RECI .NE. ' ') .AND. (RECI .NE. ',') .AND.
     :    (RECI .NE. '-'))  THEN
	i = i+1
	if (i .le. 256) go to 3
      END IF
*     ... read integer
      IEND = I-1
      ISIZE = IEND-ISTART+1
      CALL CONVRT (ISIZE,CNUM,LENTH)
      FORM = '(1I'//CNUM(1:LENTH)//')'
      READ (RECORD(ISTART:IEND),FORM,IOSTAT = IOS) LEVEL
      IF (IOS .NE. 0) THEN
          WRITE(istde,*) 'GETOLD: Unable to decode '//
     &                    RECORD(ISTART:IEND)//';'
          GOTO 1
      ENDIF
      IF (IFIRST .EQ. 0) THEN 
*       .. this is the either the first or an isolated level
         NCMIN = NCMIN+1
         IF (NCMIN .GT. NCD) THEN
             CALL RALLOC (PCCMIN,NCD,NCMIN,4)
             NCD = NCMIN
         ENDIF
         IF ((LEVEL .LT. 1) .OR. (LEVEL .GT. NCF)) THEN
            WRITE(istde,*) 'GETOLD: Serial numbers must be'
     &,              ' in the range [1,',ncf,'];'
            GOTO 1
         ENDIF
         ICCMIN(NCMIN) = LEVEL
         i = i+1
	 IF (RECI .eq. '-') ifirst = ncmin
         GO TO 2
       ELSE
*        .. the previous level was the beginning of a range
         level1 = ICCMIN(ncmin)
         number = level - level1 

         IF (number .LT. 0) THEN
            WRITE(istde,*) level1,'-',level,' not allowed'
            GOTO 1
         ENDIF

	 NCMIN = NCMIN+NUMBER
         IF (NCMIN .GT. NCD) THEN
             CALL RALLOC (PCCMIN,NCD,NCMIN,4)
             NCD = NCMIN
         ENDIF
	 DO j = 1,number
           ICCMIN(ifirst+j) =  level1 + J
	 END DO
         i = i+1
	 ifirst = 0
         GO TO 2
       END IF     

*   At least one level must be requested

    4 IF (NCMIN .EQ. 0) THEN
         WRITE(istde,*) 'GETOLD: Expecting at least one number;'
         GOTO 1
      ENDIF
*
*   Trim array to exactly the correct size
*
      IF (NCMIN .NE. NCD) THEN
         CALL RALLOC (PCCMIN,NCD,NCMIN,4)
      ENDIF
*
*   Determine NCMAX
*
      NCMAX = 0
      DO  I = 1,NCMIN
         NCMAX = MAX (NCMAX,ICCMIN(I))
      END DO
*
*   Allocate the storage for and set the weights
*
      CALL ALLOC (PWEIGH,NCMIN,8)
*
      IF (NCMIN .GT. 1) THEN
*
         IF (NDEF.NE.0) THEN
            WRITE(istde,*) 'Nonstandard level weights?'
            YES = GETYN ()
         ELSE
            YES = .FALSE.
         ENDIF
*
         IF (YES) THEN
            WRITE(istde,*) 'Assign all ASFs the same weight?'
            YES = GETYN ()
            IF (YES) THEN
               DO 5 I = 1,NCMIN
                  WEIGHT(I) = -2.0D 00
    5          CONTINUE
            ELSE
    6          WRITE(istde,*) 'Enter the (relative) weights of the'
     &,                        ncmin,' levels :'
               READ (*,*) (WEIGHT(I),I = 1,NCMIN)
               SUM = 0.0D 00
               DO 7 I = 1,NCMIN
                  IF (WEIGHT(I) .LE. 0.0D 00) THEN
                     WRITE(istde,*) 'GETOLD: Weights must exceed 0;'
                     GOTO 6
                  ELSE
                     SUM = SUM+WEIGHT(I)
                  ENDIF
    7          CONTINUE
               SUM = 1.0D 00/SUM
               DO 8 I = 1,NCF
                  WEIGHT(I) = SUM*WEIGHT(I)
    8          CONTINUE
            ENDIF
         ELSE
            SUM = 0.0D 00
            DO 9 I = 1,NCMIN
               WEIGHT(I) = -1.0D 00
    9       CONTINUE
         ENDIF
*
      ELSE
*
*   OL calculation
*
         WEIGHT(1) = -1.0D 00
*
      ENDIF
*
*   Eigenvector damping
*
      CALL ALLOC (PCDAMP,NCMIN,8)
*
      DO 10 I = 1,NCMIN
         CDAMP(I) = 0.0D 00
   10 CONTINUE
*
*   Are there any correlation functions?
*
       DO 33 I = 1,NW
         LFIX(I) = .FALSE.
   33 CONTINUE 

      WRITE(istde,*) ' Radial functions'
      CALL PRTRSL
      DO 34 I = 1,NW
         LFIX(I) = .TRUE.
   34 CONTINUE
      WRITE(istde,*) 'Enter orbitals to be varied (Updating order)'
      CALL GETRSL (INDEX,NSUBS)
      
      DO 35 I = 1,NSUBS
         LFIX(INDEX(I)) = .FALSE.
!XHH      give a big value, rather than zero to scnsty()
         scnsty(INDEX(I)) = 1.D20
   35 CONTINUE
      NFIX = NW - nsubs
      IF (NFIX .EQ. NW) THEN
*
         WRITE(istde,*) 'All subshell radial wavefunctions are fixed;'
     &,                 ' perform CI calculations with RCI92.'
!XHH
      ELSE
!         Determine orbital updating order
         NORDER = 0
         DO  I = 1,NW
            IF (.NOT. LFIX(I)) THEN
               NORDER = NORDER+1
               IORDER(I) = INDEX(NORDER)
            ENDIF
         ENDDO
      ENDIF

!XHH added a array to store the index of the correlation functions

      DO i = 1, nw
         lcorre(i) = .TRUE.
      ENDDO

      WRITE(istde,*) 'Which of these are spectroscopic orbitals?'
      CALL GETRSL (INDEX,NSUBS)
      IF (NSUBS .GT. 0) THEN
         DO 11 I = 1,NSUBS
            LOC = INDEX(I)
            IF (.NOT. LFIX(LOC)) THEN
               METHOD(LOC) = 1
               NOINVT(LOC) = .FALSE.
               ODAMP(LOC) = 0.0D 00
               lcorre(Loc) = .FALSE.
            ENDIF
   11    CONTINUE
      ENDIF

! Set NSIC. It will be non-zero if all orbitals to be varied are 
! spectroscopic orbitals

      NSIC = 4+(NW-NFIX)/4
      DO i = 1, nw
         IF ((.NOT. LFIX(i)) .AND. lcorre(i)) THEN
            NSIC = 0
            EXIT
         ENDIF
      ENDDO
*
      NSCF = 24
      NSOLV = 3
      ORTHST = .TRUE.
*
*   Make the allocation for the auxiliary vector required
*   by SUBROUTINE NEWCO
*
      CALL ALLOC (PNTRWT,NCMIN,8)
*
      RETURN
      END
