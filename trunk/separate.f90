subroutine separate()
use sfvars, only: natom, atpos, attyp, atsel, nimage, xlo,xhi,ylo,yhi,zlo,zhi, &
nsel, list, istr, iend
implicit none
   !----------------------------------------------------------------------------
   integer             :: i, j, iou, ii
   character (len=100) :: fname
   !----------------------------------------------------------------------------
   write(*,'(/,10x,"Please input the filename to output the result: ",$)')
   read(*, *) fname
   iou = 10
   open(iou, file=fname, status='unknown', action='write')
   do ii = 1, nsel
      i = list(ii)
      call header()
      do j = istr, iend
         write(iou, 100) j, attyp(i), atpos(:,i,j)
      enddo
   enddo
   close(iou)
100 format(I6,1x,I2,3(1x,F20.10))
   !----------------------------------------------------------------------------
contains
   !----------------------------------------------------------------------------
   subroutine header( )
   implicit none
      !-------------------------------------------------------------------------
      write(iou, '("ITEM: TIMESTEP")')
      write(iou, '(I10)') i
      write(iou, '("ITEM: NUMBER OF ATOMS")')
      write(iou, '(I10)') iend-istr+1
      write(iou, '("ITEM: BOX BOUNDS")')
      write(iou, '(F15.6,1x,F15.6)') xlo, xhi
      write(iou, '(F15.6,1x,F15.6)') ylo, yhi
      write(iou, '(F15.6,1x,F15.6)') zlo, zhi
      write(iou, '("ITEM: ATOMS")')
      !-------------------------------------------------------------------------
   return
   end subroutine
   !----------------------------------------------------------------------------
end subroutine
      
