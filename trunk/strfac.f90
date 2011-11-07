! Subroutine to calculate the structure factor for selected atoms of one frame
! Note that atpos is in fractional coordinate
subroutine strfaccal
use sfvars
implicit none
   !----------------------------------------------------------------------------
   integer :: l, m, n, i, II
   real(q) :: kx, ky, kz
   complex            :: ctmp, jr(3)
   complex, parameter :: img = (0.D0,1.D0)*tpi
   !----------------------------------------------------------------------------
   !
   do II = 1, nsel
      i  = list(II)
      jr = OneImg(:,i)*img
      do l = 1, Nx
         kx = dble(l-1)
         qxrx(II,l) = exp( kx*jr(1) )
      enddo
      do m = 1, Ny
         ky = dble(m-1)
         qyry(II,m) = exp( ky*jr(2) )
      enddo
      do n = 1, Nz
         kz = dble(n-1)
         qzrz(II,n) = exp( kz*jr(3) )
      enddo
   enddo
   !
   do n = 1, Nz
   do m = 1, Ny
   do l = 1, Nx
      ctmp = (0.D0,0.D0)
      do II = 1, nsel
         ctmp = ctmp + qxrx(II,l)*qyry(II,m)*qzrz(II,n)
      enddo
      sqcur(l,m,n) = ctmp
   enddo
   enddo
   enddo
   !
   sqcur = sqcur / sqrt(dble(nsel))
   !----------------------------------------------------------------------------
return
end subroutine
