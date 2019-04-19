!
! Sometimes the 'mcp.{kmax}' file has zero record due to the fact
! that its contents are not flushed out when mcp crashes (for
! example, file system full...). This patch file helps restarting
! mcp.
! Before running this utility, the "magic numbers" should be 
! obtained from a restart run of mcp (they are from lodres.f).
!
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL DIAG,LFRDRT
      DIAG = .FALSE.
      LFRDRT = .FALSE.
      WRITE (41) 'MCP','UNSORTED'
      WRITE (41) 103, 6839, 44
      WRITE (41) DIAG, 0, LFRDRT
      END
