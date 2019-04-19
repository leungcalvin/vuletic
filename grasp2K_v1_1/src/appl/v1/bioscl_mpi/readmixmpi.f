************************************************************************
*                                                                      *
      SUBROUTINE READMIX(NAME,INPCI)
*                                                                      *
*   Open and read the mixing coefficent files                          *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

Cww      INTEGER PNTRIQ, PNTIQR
      POINTER(PNTRIQ,RIQDUMMY)
      POINTER(PNTIQR,IQRDUMMY)

      CHARACTER*128 NAME(2)
      CHARACTER*6 G92MIX

      POINTER(PNEVALII,EVALII(1))
      POINTER(PNEVECII,EVECII(1))
      POINTER(PNIVECII,IVECII(1))
      POINTER(PIATJPII,IATJPOII(1)),(PIASPAII,IASPARII(1))

      POINTER(PNEVALFF,EVALFF(1))
      POINTER(PNEVECFF,EVECFF(1))
      POINTER(PNIVECFF,IVECFF(1))
      POINTER(PIATJPFF,IATJPOFF(1)),(PIASPAFF,IASPARFF(1))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF1R/NELECR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB2R/NCFR,NWR,PNTIQR

      COMMON/DEF1II/EMNII,IONCTYII,NELECII,ZII
     :      /EIGVALII/EAVII,PNEVALII
     :      /EIGVECII/PNEVECII
     :      /ORB2II/NCFII,NWII
     :      /PRNTII/NVECII,PNIVECII,NVECMXII
     :      /SYMAII/PIATJPII,PIASPAII

      COMMON/DEF1FF/EMNFF,IONCTYFF,NELECFF,ZFF
     :      /EIGVALFF/EAVFF,PNEVALFF
     :      /EIGVECFF/PNEVECFF
     :      /ORB2FF/NCFFF,NWFF
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF
     :      /SYMAFF/PIATJPFF,PIASPAFF
*
*   Read the initial state mixing file    
*
      J = INDEX(NAME(1),' ')
      IF (INPCI.EQ.0) THEN
        OPEN (UNIT = 68,FILE=NAME(1)(1:J-1)//'.cbm',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ELSE
        OPEN (UNIT = 68,FILE=NAME(1)(1:J-1)//'.bm',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ENDIF
      READ(68,IOSTAT=IOS) G92MIX
      IF ((IOS.NE.0).OR.(G92MIX.NE.'G92MIX')) THEN
        WRITE(*,*) 'Not a GRASP mixing file'
        STOP
      ENDIF

      READ(68) N11,N12,N13
      IF ((N11.NE.NELECII).OR.(N12.NE.NCFII).OR.(N13.NE.NWII)) THEN
        PRINT *, 'This MIXing Coefficients File is not'
        PRINT *, 'appropriate for the initial state'
        STOP
      ENDIF
*
      READ(68) NVECII
      CALL ALLOC(PNEVALII,NVECII,8)
      CALL ALLOC(PNEVECII,NCFII*NVECII,8)
      CALL ALLOC(PNIVECII,NVECII,4)
      CALL ALLOC(PIATJPII,NVECII,4)
      CALL ALLOC(PIASPAII,NVECII,4)
      READ(68) (IVECII(I),I=1,NVECII)
      READ(68) (IATJPOII(I),IASPARII(I),I=1,NVECII)
      READ(68) EAVII,(EVALII(I),I=1,NVECII)

C      DO 10 I=1,NVECII
C         EVALII(I)=EAVII+EVALII(I)
C10    CONTINUE
      READ(68) ((EVECII(I+(J-1)*NCFII),I=1,NCFII),J=1,NVECII)
*
*   Close the initial state mixing  file
*
      CLOSE(68)
*
*   Read the final state mixing file    
*
      J = INDEX(NAME(2),' ')
      IF (INPCI.EQ.0) THEN
        OPEN (UNIT = 78,FILE=NAME(2)(1:J-1)//'.cbm',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ELSE 
        OPEN (UNIT = 78,FILE=NAME(2)(1:J-1)//'.bm',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ENDIF

      READ(78,IOSTAT=IOS) G92MIX
      IF ((IOS.NE.0).OR.(G92MIX.NE.'G92MIX')) THEN
        WRITE(*,*) 'Not a GRASP mixing file'
        STOP
      ENDIF

      READ(78) N11,N12,N13
      IF ((N11.NE.NELECFF).OR.(N12.NE.NCFFF).OR.(N13.NE.NWFF)) THEn
        PRINT *, 'This MIXing Coefficients File is not'
        PRINT *, 'appropriate for the final state'
        STOP
      ENDIF
*
      READ(78) NVECFF
      CALL ALLOC(PNEVALFF,NVECFF,8)
      CALL ALLOC(PNEVECFF,NCFFF*NVECFF,8)
      CALL ALLOC(PNIVECFF,NVECFF,4)
      CALL ALLOC(PIATJPFF,NVECFF,4)
      CALL ALLOC(PIASPAFF,NVECFF,4)
      READ(78) (IVECFF(I),I=1,NVECFF)
      READ(78) (IATJPOFF(I),IASPARFF(I),I=1,NVECFF)
      READ(78) EAVFF,(EVALFF(I),I=1,NVECFF)

C      DO 20 I=1,NVECFF
C         EVALFF(I)=EAVFF+EVALFF(I)
C20    CONTINUE
      READ(78) ((EVECFF(I+(J-1)*NCFFF),I=1,NCFFF),J=1,NVECFF)
*
*   Close the initial state mixing  file
*
      CLOSE(78)

      RETURN
      ENd
