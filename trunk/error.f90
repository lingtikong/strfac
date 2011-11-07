!***********************************************************************
!***            Subroutine for warning and error                     ***
!***********************************************************************
module errmesg
implicit none
   !----------------------------------------------------------------------------
   character(len=10)  :: subname
   character(len=100) :: errinfo
   !
contains
   subroutine warn(ierr)
   implicit none
      !--------------------------------------------------------------------
      integer :: ierr
      !--------------------------------------------------------------------
      write(*, '(/, 10X, 60("W") )')
      write(*, '( 10X, "WWW", 14X, "WARNING from ", A12, 15X, "WWW")') subname
      write(*, '( 10X, "WWW", 54X, "WWW")')
      write(*, '( 10X, "WWW", 2X, A50, 2X, "WWW" )'   ) trim( errinfo(1:50) )
      if ( trim( errinfo(51:) ).ne.'' ) &
      &  write(*, '( 10X, "WWW", 2X, A50, 2X, "WWW" )'  ) trim( errinfo(51:) )
      write(*, '( 10X, "WWW", 2X, "Error Code: ", I4, 36X, "WWW" )' ) ierr
      write(*, '( 10X, "WWW", 54X, "WWW", /, 10x, 60("W") )')
      !--------------------------------------------------------------------
   return
   end subroutine warn
   !
   subroutine error(ierr)
   implicit none
      !--------------------------------------------------------------------
      integer :: ierr
      !-------------------------------------------------------------------
      if ( ierr.eq.0 ) return
      !
      write(*, '(/, 10X, 60("E") )')
      !
      write(*, '( 10X, "EEE", 16X, "ERROR from ", A12, 15X, "EEE")') subname
      write(*, '( 10X, "EEE", 54X, "EEE")')
      write(*, '( 10X, "EEE", 2X, A50, 2X, "EEE" )'   ) trim( errinfo(1:50) )
      if ( trim( errinfo(51:) ).ne.'' ) &
      &  write(*, '( 10X, "EEE", 2X, A50, 2X, "EEE" )' ) trim( errinfo(51:) )
      write(*, '( 10X, "EEE", 2X, "Error Code: ", I4, 36X, "EEE" )' ) ierr
      write(*, '( 10X, "EEE", 54X, "EEE", /, 10X, 60("E") )')
      !
      if ( ierr.gt.0 ) stop
      !---------------------------------------------------------------------
   return
   end subroutine error
end module
