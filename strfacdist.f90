! Module to evaluate the static structure factor distribution at a given k
! Note that atpos is in fractional coordinate
module SFDist
use sfvars
use errmesg
use iounits
implicit none
   private
   public :: StrFacDist
   !----------------------------------------------------------------------------
   real(q) :: qk(3), rprd(3)
   complex, allocatable :: Sq(:,:), dist(:)
   integer, allocatable :: nhit(:)
   !----------------------------------------------------------------------------
contains
   !
subroutine StrFacDist()
implicit none
   !----------------------------------------------------------------------------
   integer :: Iimg, nimg, i, II, Im
   complex :: ctmp, qi(3)
   complex, parameter :: img = (0.D0,1.D0)*tpi
   !
   real(q) :: ri, dr, rdr, rmin, rmax, Rc(3), RiC(3)
   integer :: npt, disttype, idx
   !----------------------------------------------------------------------------
   write(*,'(/,10x,"Please input the vector q in units of (2pi/Lx 2pi/Ly 2pi/Lz): ",$)')
   read(*,*,iostat=ioerr) qk
   if ( ioerr.ne.0 ) then
      write(errinfo, '("Wrong input!")')
      call warn(ioerr)
      return
   endif
   !
   nimg = (iend-istr)/inc+1
   if (allocated(Sq)) deallocate(Sq)
   allocate(Sq(nsel,nimg))
   !
   qi = img * qk * 8.D0*atan(1.D0)
   !
   Im = 1
   do Iimg = istr, iend, inc
      OneImg = atpos(:,:,Iimg)
      !
      do II = 1, nsel
         i  = list(II)
         ctmp = OneImg(1,i)*qi(1) + OneImg(2,i)*qi(2) + OneImg(3,i)*qi(3)
         ctmp = exp(ctmp)
         Sq(II, Im) = ctmp
      enddo
      Im = Im + 1
   enddo
   !
   disttype = 3
   write(*, '(/,10x,"***result might be unreliable for non-orthogonal box***")')
   write(*, '(10x,"Please select the distribution type")')
   write(*, '(12x,"1. <S^2(k)> along x;")')
   write(*, '(12x,"2. <S^2(k)> along y;")')
   write(*, '(12x,"3. <S^2(k)> along z;")')
   write(*, '(12x,"4. radial distribution of <S^2(k)>;")')
   write(*, '(10x,"Your choice [3]: ")')
   read(*,'(A)', iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') then
      read(oneline,*,iostat=ioerr) disttype
      if (ioerr.ne.0) disttype = 3
   endif
   write(*, '(10x,"Your selected:",I2)') disttype
   npt = 100
   write(*, '(10x,"Please input the # of points along selected direction [100]:", $)')
   read(*,'(A)', iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') then
      read(oneline,*,iostat=ioerr) npt
      if (ioerr.ne.0) npt = 3
   endif
   if (npt < 1) return
   !
   if (allocated(dist)) deallocate(dist)
   if (allocated(nhit)) deallocate(nhit)
   allocate(nhit(npt+1), dist(npt+1))
   nhit = 0
   dist = (0.D0,0.D0)
   !
   select case (disttype)
   case ( :3 )
      rmin = 0.D0
      rmax = 1.D0
      dr   = (rmax-rmin)/dble(npt)
      rdr  = 1.D0/dr

      Im = 1
      do Iimg = istr, iend, inc
         OneImg = atpos(:,:,Iimg)
         do II = 1, nsel
            i  = list(II)
            ri = OneImg(disttype,i)
            ri = ri - floor(ri)
            idx= int(ri*rdr) + 1
            nhit(idx) = nhit(idx) + 1
            dist(idx) = dist(idx) + Sq(II,Im)
         enddo
         Im = Im + 1
      enddo
      !
      prd  = box(:,istr)
      rmax = prd(disttype)
      dr   = (rmax-rmin)/dble(npt)
      rdr  = 1.D0/dr
      !
   case ( 4 )
      !
      write(*,'(/,10x,"***All distance in this part should be in cartesian***")')
      write(*,'(10x,"Please input the coordinate of the origin [0 0 0]: ", $)')
      read(*, '(A)', iostat=ioerr) oneline
      if (ioerr.ne.0) return
      if (oneline.ne.'') then
         read(oneline, *, iostat=ioerr) Rc
         if (ioerr.ne.0) return
      endif
      prd  = box(:,istr)
      rmin = 0.D0
      rmax = 0.5D0*minval(prd(1:3))
      write(*,'(10x,"Please input the maximum distance to get dist [",F5.2,"]: ", $)') rmax
      read(*, '(A)', iostat=ioerr) oneline
      if (ioerr.ne.0) return
      if (oneline.ne.'') then
         read(oneline, *, iostat=ioerr) rmax
         if (ioerr.ne.0) return
      endif
      dr   = (rmax-rmin)/dble(npt)
      rdr  = 1.D0/dr

      Im = 1
      do Iimg = istr, iend, inc
         OneImg = atpos(:,:,Iimg)
         prd    = box(:,Iimg)
         rprd   = 1.D0/prd(1:3)
         do II = 1, nsel
            i  = list(II)
            RiC = OneImg(:,i) - Rc
            RiC = Ric - nint(Ric*rprd)*prd(1:3)
            ri = sqrt(sum(RiC*RiC))
            idx= int(ri*rdr) + 1
            nhit(idx) = nhit(idx) + 1
            dist(idx) = dist(idx) + Sq(II,Im)
         enddo
         Im = Im + 1
      enddo
      !
   case default
      return
   end select
   !
   where (nhit>0) dist = dist/dble(nhit)
   dist = dist*conjg(dist)
   !
   ! write out the result
   outfile = "SqDist.dat"
   write(*,'(10x,"Please input the filename to output S(q) distribution [",A,"]: ",$)') trim(outfile)
   read(*,'(A)', iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') then
      read(oneline,*,iostat=ioerr) outfile
      if (ioerr.ne.0) outfile = "SqDist.dat"
   endif
   open(iotmp, file=outfile, action='write')
   write(iotmp,'("# |<S(q)>|^2 at q = [",3F8.4,"].")') qk
   write(iotmp,'("# r   |<S(q)>|^2")')
   ri = 0.5D0*dr
   do i = 1, npt
      write(iotmp,'(F12.6,1x,F20.10)') ri, real(dist(i))
      ri = ri + dr
   enddo
   close(iotmp)
   !
   if (allocated(Sq))   deallocate(Sq)
   if (allocated(dist)) deallocate(dist)
   if (allocated(nhit)) deallocate(nhit)
   !----------------------------------------------------------------------------
return
end subroutine
   !----------------------------------------------------------------------------
end module
