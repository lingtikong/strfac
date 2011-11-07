module sfvars
use prec
implicit none
   !----------------------------------------------------------------------------
   integer, parameter :: max_type = 50
   !----------------------------------------------------------------------------
   integer :: ioerr, istr, iend, inc = 1, ngr, job = 1, NimgMax
   integer :: Nx, Ny, Nz, natom, nsel, ntype, image, nimage, nimgused
   !----------------------------------------------------------------------------
   real(q) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
   real(q) :: grbase, grwt
   real(q) :: Lx, Ly, Lz, prd(6), qvec(3), rmax, dr, rdr, halfdr, vol
   integer :: plane(3) = 1, typsrc, typdes, nsrc, ndes, ir2r
   real(q), allocatable :: atpos(:,:,:), OneImg(:,:), box(:,:)
   integer, allocatable :: attyp(:), atsel(:), list(:), timestep(:), Indice(:)
   integer, allocatable :: srclist(:), deslist(:), nncur(:)
   integer, allocatable :: wrap(:,:,:), wrapone(:,:)
   logical              :: lwrapinfo = .false.
   !
   real(q), allocatable :: sqsum(:,:,:), grsum(:), nnsum(:)
   complex, allocatable :: sqcur(:,:,:), qxrx(:,:), qyry(:,:), qzrz(:,:)
   !----------------------------------------------------------------------------
   character(len=256)   :: oneline
   character(len=256)   :: infile, outfile
   !----------------------------------------------------------------------------
   real(q), parameter   :: tpi = 2.D0*3.1415926536D0
   !----------------------------------------------------------------------------
   real(q)              :: tmbeg, tmend, tmread, tmcal
   !----------------------------------------------------------------------------
   !
contains
   subroutine free_all
   implicit none
      if (allocated(atpos   )) deallocate( atpos   )
      if (allocated(OneImg  )) deallocate( OneImg  )
      if (allocated(box     )) deallocate( box     )
      if (allocated(attyp   )) deallocate( attyp   )
      if (allocated(atsel   )) deallocate( atsel   )
      if (allocated(list    )) deallocate( list    )
      if (allocated(sqsum   )) deallocate( sqsum   )
      if (allocated(grsum   )) deallocate( grsum   )
      if (allocated(nnsum   )) deallocate( nnsum   )
      if (allocated(nncur   )) deallocate( nncur   )
      if (allocated(sqcur   )) deallocate( sqcur   )
      if (allocated(qxrx    )) deallocate( qxrx    )
      if (allocated(qyry    )) deallocate( qyry    )
      if (allocated(qzrz    )) deallocate( qzrz    )
      if (allocated(timestep)) deallocate( timestep)
   end subroutine
end module
