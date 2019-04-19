************************************************************************
*                                                                      *
      SUBROUTINE LODMIX
*                                                                      *
*   Determines the eigenpairs required;  this information is written   *
*   to the head of the  .mix  file.                                    *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, RALLOC.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 RECORD
      CHARACTER*7 FORM
      CHARACTER*3 CNUM
      CHARACTER*1 RECI
*
      POINTER (PNIVEC,IVEC(1))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
      COMMON/iounit/istdi,istdo,istde
*
*   Make the initial allocation for the level list vector
*
      NVECD = 1
      CALL ALLOC (PNIVEC,NVECD,4)
*
*   Entry message
*
    1 WRITE(istde,*) 'Enter the serial number(s) of'
     &, ' the level(s) to be calculated:'
*
*   Read and parse the list of levels
*
      READ (*,'(A)') RECORD
*
*   Initialise NVEC
*
      NVEC = 0
*
*   Parse RECORD from left to right
*
      IFIRST = 0
      ISTART = 1
      I = 1
*
*   .. skip the blanks and commas (this implementation allows input to 
*      start with blanks
    2 RECI = RECORD(I:I)

      IF ((RECI .NE. ' ') .AND. (RECI .NE. ',')) THEN

!       skip all other characters except 1-9
!      IF ((RECI .GE. '1') .AND. (RECI .LE. '9')) THEN
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
     :    (RECI .NE. '-')) THEN
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
         WRITE(istde,*) 'LODMIX: Unable to decode '
     :                //RECORD(ISTART:IEND)//';'
          GOTO 1
      ENDIF
      IF (IFIRST .EQ. 0) THEN 
*       .. this is the either the first or an isolated level
         NVEC = NVEC+1
         IF (NVEC .GT. NVECD) THEN
             CALL RALLOC (PNIVEC,NVECD,NVEC,4)
             NVECD = NVEC
         ENDIF
         IF ((LEVEL .LT. 1) .OR. (LEVEL .GT. NCF)) THEN
            WRITE(istde,*) 'LODMIX: Serial numbers must be'
     &              , ' in the range [1,',NCF,'];'
            GOTO 1
         ENDIF
         IVEC(NVEC) = LEVEL
         i = i+1
	 IF (RECI .eq. '-') ifirst = NVEC
         GO TO 2
      ELSE
*        .. the previous level was the beginning of a range
         level1 = IVEC(ifirst)
         number = level - level1
         IF (number .LT. 0) THEN
            WRITE(istde,*) level1,'-',level,' not allowed'
            GOTO 1
         ENDIF
	 NVEC = NVEC+number
         IF (NVEC .GT. NVECD) THEN
             CALL RALLOC (PNIVEC,NVECD,NVEC,4)
             NVECD = NVEC
         ENDIF
	 DO j = 1,number
           IVEC(ifirst+j) =  level1 + J
	 END DO
         i = i+1
	 ifirst = 0
         GO TO 2
      END IF     
*
*   At least one level must be requested
*
    4 IF (NVEC .EQ. 0) THEN
         WRITE(istde,*) 'LODMIX: Expecting at least one number;'
         GOTO 1
      ENDIF
*      print *, (ivec(j),j=1,nvec)
*      print *, 'nvec, nvecd', nvec, nvecd
*      goto 1
*   Trim array to exactly the correct size
*
      IF (NVEC .NE. NVECD) THEN
         CALL RALLOC (PNIVEC,NVECD,NVEC,4)
      ENDIF
*
*   Determine NVECMX
*
      NVECMX = 0
      DO 5 I = 1,NVEC
         NVECMX = MAX (NVECMX,IVEC(I))
    5 CONTINUE
*
*   Write header data to the  .mix  file
*
      WRITE (25) NELEC,NCF,NW
      WRITE (25) NVEC
      WRITE (25) (IVEC(I),I = 1,NVEC)
*
      RETURN
      END
************************************************************************
