************************************************************************
*                                                                      *
      SUBROUTINE FNDBLK (LPRINT)
*                                                                      *
*   The CSFs are examined to determine the block structure or the
*   matrix. Information is  returned in COMMON/HBLOCK/   *
*   NBLOCK  - number of blocks                                         *
*   NELBLK(.) - number of CSFs in each block
*   IBLOCK(.) - index of block for each CSF
*
*   Call(s) to: ALLOC, DALLOC, RALLOC.                                 *
*                                                                      *
*   Written by Charlotte F. Fischer
*                                         Last revision: 13 May 1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNWP = 30)
      LOGICAL LPRINT
      CHARACTER*256 record
      CHARACTER*8 cnum
*
      POINTER (PIBLOC,IBLOCK(1))
      POINTER (PNELBL,NELBLK(1))
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNTJQS,JQSDUMMY)
*
      INTEGER*4 JCUPA
      POINTER (PNJCUP,JCUPA(NNNWP,1))
*
      COMMON/HBLOCK/NBLOCK,PIBLOC,PNELBL
     :      /ORB2/NCF,NW,PNTRIQ
     :      /STAT/PNTJQS,PNJCUP

      POINTER (PIDBLK,IDBLK(1)), (PNEL, NEL(1))
*
*   Entry message
*
      WRITE (*,*) 'Calling FNDBLK...'
*
*   Allocate storage for IBLOCK, NELBL and initialize
*
      CALL ALLOC (PIBLOC,NCF,4)
      NELD = 10 
      CALL ALLOC (PIDBLK, NELD, 4)
      CALL ALLOC (PNEL, NELD, 4)
      NBLOCK = 0
      DO 1 I = 1,NCF
         IBLOCK(I) = 0
    1 CONTINUE
*
*     .. classify each CSF
*
      DO IC = 1,NCF
*    
*	.. find J and parity
	jp = ibits(jcupa((NNNW-1)/4+1,IC),8*MOD(NNNW-1,4),8)
*	  .. is this block already identified
	  ipos = 0
	  DO i = 1,nblock
	    IF (jp .eq. idblk(i)) THEN
	      ipos = i
	      EXIT
	    END IF
	  END DO
	  IF (ipos .ne. 0) THEN
*           .. J and parity already defined
	    iblock(IC) = ipos
	    nel(ipos) = nel(ipos) + 1
	  ELSE
*           .. a new block has been found
	    nblock = nblock + 1
	    IF (nblock .gt. NELD) THEN
*             .. list must be longer
	      newsize = neld + neld/2
	      CALL ralloc (pnel, neld, newsize, 4)
	      CALL ralloc (pidblk, neld, newsize, 4)
	      neld = newsize
	    END IF
	    idblk(nblock) = jp
	    iblock(ic) = nblock
	    nel(nblock) =1
	  END IF
	END DO
*
*       .. save nel in block of proper size
	CALL alloc (pnelbl, nblock, 4)
	DO i = 1,nblock
	  nelblk(i) = nel(i)
	END DO
	CALL dalloc (PNEL)
	CALL dalloc (PIDBLK)
*
*   Print summary
*
      IF (NBLOCK .GT. 1) THEN
         WRITE (*,301) NBLOCK
         WRITE (*,302) (NELBLK(I),I = 1,NBLOCK)
         IF (LPRINT) THEN
            WRITE (99,303) (IBLOCK(I),I = 1,NCF)
         ENDIF
      ENDIF
*
      RETURN
*
  301 FORMAT (' Matrix is reducible; number of blocks is ',1I6,'.')
  302 FORMAT (1X,' Number of elements in each block:'/1X,13(1X,1I6))
  303 FORMAT (1X,' Index array:' /44(1X,1I2))
*
      END

