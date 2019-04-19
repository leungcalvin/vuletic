************************************************************************
*                                                                      *
      SUBROUTINE PRSRSL (NFILE,ID)
*                                                                      *
*   READs and parses a list of subshell labels on unit NFILE to load   *
          include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
*   COMMON  blocks  /ORB4/,  /ORB5/, and /ORB10/; the value of NW in   *
*   COMMON/ORB2/ is incremnted by 1 for each new subshell label; the   *
*   labels are delimited either by blanks or commas.                   *
*                                                                      *
*   Call(s) to: [LIB92] CONVRT.                                        *
*                                                                      *
*   Written by Farid A Parpia               Last revised 14 Sep 1992   *
*                                                                      *
************************************************************************
      COMMON/iounit/istdi,istdo,istde
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
CGG      CHARACTER*500 RECORD, RECORD2,REC
      CHARACTER*1000 RECORD, RECORD2,REC
      CHARACTER*6 FORM
      CHARACTER*2 CTEGER,NH,SYM,SYMLST
      CHARACTER*1 RECI
      LOGICAL NEWREC
*
      DIMENSION SYMLST(19)
*
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
*
      DATA SYMLST/'s ','p-','p ','d-','d ','f-','f ','g-','g ','h-',
     :            'h ','i-','i ','k-','k ','l-','l ','m-','m'/
*
*   Read the records
*
      NEWREC = .FALSE.
      READ (NFILE,'(A)') RECORD
      IF (ID.EQ.2) THEN
         READ (NFILE,'(A)') RECORD2
         IF (RECORD2(1:3).EQ.'CSF') THEN
            BACKSPACE (UNIT = NFILE)
         ELSE
            NEWREC = .TRUE.
         ENDIF
      ENDIF
*
*   Parse RECORD from left to right
*
      ISTART = 0
CGG
      IOS = 0
CGG
      I = 1
    1 RECI = RECORD(I:I)
      IF ((RECI .NE. ' ') .AND. (RECI .NE. ',')) THEN
         IF (ISTART .EQ. 0) ISTART = I
      ELSE
         IF (ISTART .NE. 0) THEN
            IEND = I-1
            RECI = RECORD(IEND:IEND)
            IF (RECI .EQ. '-') THEN
c              READ (RECORD(IEND-1:IEND),'(1A2)',IOSTAT = IOS) SYM
               SYM = RECORD(IEND-1:IEND)
               IF (IOS .NE. 0) THEN
                 WRITE(istde,*) 'PRSRSL: Symmetry ',RECORD(IEND-1:IEND)
     &  , ' could not be decoded.'
                  STOP
               ENDIF
               IEND = IEND-2
            ELSE
               SYM(2:2) = ' '
c              READ (RECI,'(1A1)',IOSTAT = IOS) SYM(1:1)
               SYM(1:1) = RECI
               IF (IOS .NE. 0) THEN
                  WRITE(istde,*) 'PRSRSL: Symmetry ',RECI,' could not'
     & , ' be decoded.'
                  STOP
               ENDIF
               IEND = IEND-1
            ENDIF
            DO 2 II = 1,19
               IF (SYM .EQ. SYMLST(II)) THEN
                  NW = NW+1
                  IF (NW .GT. NNNW) THEN
                     WRITE(istde,*) 'PRSRSL: Number of subshells '
     & , 'exceeds allocation: plant NW was set to',NNNW 
                     STOP
                  ENDIF
                  NH(NW) = SYM
                  IIB2 = II/2
                  NKL(NW) = IIB2
                  IF (MOD (II,2) .EQ. 1) THEN
                     NAK(NW) = -IIB2-1
                     NKJ(NW) = II
                  ELSE
                     NAK(NW) =  IIB2
                     NKJ(NW) = II-1
                  ENDIF
                  CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
                  FORM = '(1I'//CTEGER(1:LTEGER)//')'
                  READ (RECORD(ISTART:IEND),FMT = FORM,IOSTAT = IOS)
     :               NP(NW)
                  IF (IOS .NE. 0) THEN
                     WRITE(istde,*) 'PRSRSL: Principal quantum number ',
     :                        RECORD(ISTART:IEND)
     &, ' could not be decoded.'
                     STOP
                  ENDIF
                  GOTO 3
               ENDIF
    2       CONTINUE
            WRITE(istde,*) 'PRSRSL: Symmetry ',SYM,
     &                     ' could not be decoded.'
            STOP
    3       ISTART = 0
         ENDIF
      ENDIF
*
CGG      IF (I .LT. 500) THEN
      IF (I .LT. 1000) THEN
         I = I+1
         GOTO 1
      ENDIF
      IF (NEWREC) THEN 
         RECORD = RECORD2
         NEWREC = .FALSE.
         I = 1
         ISTART = 0
         GOTO 1
      ENDIF
*
      RETURN
      END
