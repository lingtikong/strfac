subroutine SISF()
use sfvars
implicit none
   !----------------------------------------------------------------------------
   real(q) :: qk, Fcur
   real(q) :: Rij(3), Rjk(3), Rik(3), r2
   integer :: i, j, ii, k, it, jt, is
   real(q), allocatable :: Fqt(:)
   integer, allocatable :: stepdiff(:)
   character (len=100)  :: fout
   character(len=1)     :: signals(3) = (/ '\', '|', '/' /)
   !----------------------------------------------------------------------------
   write(*,'(/,10x,"Please input the value for q: ", $)')
   read(*, *, iostat=ioerr) qk
   if (ioerr.ne.0) return
   !
   if (allocated(Fqt)) deallocate(Fqt)
   if (allocated(stepdiff)) deallocate(stepdiff)
   allocate( Fqt(nimgused), stepdiff(nimgused) )
   Fqt    = 0.D0
   !
   write(*, '(10x,"Computing ...", $)')
   it = 0
   is = 0
   do i = istr, iend, inc
      it = it + 1
      jt = it + 1
      do j = i+inc, iend, inc
         jt = jt + 1
         !----------------------------------------------------------------------
         Fcur = 0.D0
         do ii = 1, nsel
            k  = list(ii)
            Rjk = atpos(:,k,j) + dble(wrap(:,k,j))
            Rjk(1) = Rjk(1) * box(1,j) + Rjk(2)*box(4,j) + Rjk(3)*box(5,j)
            Rjk(2) = Rjk(2) * box(2,j) + Rjk(3)*box(6,j)
            Rjk(3) = Rjk(3) * box(3,j)
            Rik = atpos(:,k,i) + dble(wrap(:,k,i))
            Rik(1) = Rik(1) * box(1,i) + Rik(2)*box(4,i) + Rik(3)*box(5,i)
            Rik(2) = Rik(2) * box(2,i) + Rik(3)*box(6,i)
            Rik(3) = Rik(3) * box(3,i)
            RIJ = Rjk - Rik
            r2  = sum(RIJ*RIJ)
            Fcur = Fcur + cos(qk*sqrt(r2))
         enddo
         Fqt(jt-it) = Fqt(jt-it) + Fcur /dble(nsel)
         is = mod(is,3) + 1
         write(*,'(2A1,$)') char(8), signals(is)
      enddo
      stepdiff(it) = timestep(i) - timestep(istr)
   enddo
   write(*,'("Done!")')
   do i = 1, nimgused
      Fqt(i) = Fqt(i)/dble(nimgused-i+1)
   enddo
   Fqt(1) = 1.D0
   !
   write(*,'(10x,"Please input the output file name [sisf.dat]: ", $)')
   read(*,'(A)', iostat=ioerr) fout
   if (ioerr.ne.0.or.fout.eq."") fout = "sisf.dat"
   !
   open(10, file=fout, action='write')
   write(10,'("#stepdiff  SISF")')
   do i = 1, nimgused
      write(10,'(I10,2x,F20.10)') stepdiff(i), Fqt(i)
   enddo
   close(10)
   !----------------------------------------------------------------------------
end subroutine
