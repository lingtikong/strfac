! Subroutine to calculate the pair correlation for selected atoms of one frame
! Note that atpos is in fractional coordinate
subroutine paircorr
use sfvars
implicit none
   !----------------------------------------------------------------------------
   integer :: II, JJ, i, j, ip
   real(q) :: Rij(3), r, posi(3)
   !----------------------------------------------------------------------------
   nncur = 0
   do II = 1, nsel-1
      i    = list(II)
      posi = OneImg(:,i)
      do JJ = II+1, nsel
         j = list(JJ)
         !
         Rij = OneImg(:,j) - posi
         Rij = (Rij - NINT(Rij))
         Rij(1) = Rij(1) * prd(1) + Rij(2)*prd(4) + Rij(3)*prd(5)
         Rij(2) = Rij(2) * prd(2) + Rij(3)*prd(6)
         Rij(3) = Rij(3) * prd(3)
         r   = sqrt(sum(Rij*Rij))
         !
         ip = int(r*rdr) + 1
         if (ip.le.ngr) nncur(ip) = nncur(ip) + 1
         !
      enddo
   enddo
   !----------------------------------------------------------------------------
return
end subroutine
