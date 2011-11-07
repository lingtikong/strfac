! Subroutine to calculate the structure factor for selected atoms of one frame
! Note that atpos is in fractional coordinate
subroutine strfacOne
use sfvars
implicit none
   !----------------------------------------------------------------------------
   integer            :: i, II
   complex            :: jr(3)
   complex, parameter :: img = (0.D0,1.D0)*tpi
   !----------------------------------------------------------------------------
   !
   do II = 1, nsel
      i  = list(II)
      jr = OneImg(:,i)*img
      sqcur(image,1,1) = sqcur(image,1,1) + exp(sum(jr*qvec))
   enddo
   !----------------------------------------------------------------------------
return
end subroutine
