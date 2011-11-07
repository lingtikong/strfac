! Subroutine to output structure factor
subroutine outputsf
use sfvars
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer, parameter :: npt = 120
   !----------------------------------------------------------------------------
   integer :: i,j,k, iq, nsq(npt)
   real(q) :: qr, dq, rdq, sq(npt)
   !----------------------------------------------------------------------------
   outfile = 'strfac.dat'
   write(*,'(/,10x,"Please input the output file name for raw data   [strfac.dat]: ",$)')
   read(*,'(A)',iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') outfile = trim(oneline)
   !
   nsq = 0
   sq  = 0.D0
   qvec = (/ dble(Nx)/Lx, dble(Ny)/Ly, dble(Nz)/Lz /) * tpi
   qr   = sqrt( sum(qvec*qvec) )
   dq   = qr / dble(npt)
   rdq  = 1.D0/dq
   !
   open(ioout, file=outfile, status='unknown', action='write')
   write(ioout, 100)
   do k = 1, Nz
   do j = 1, Ny
   do i = 1, Nx
      qvec = (/ dble(i-1)/Lx, dble(j-1)/Ly, dble(k-1)/Lz /) *tpi
      qr   = sqrt( sum(qvec*qvec) )
      iq   = int(qr*rdq)+1
      nsq(iq) = nsq(iq) + 1
      sq(iq)  = sq(iq)  + sqsum(i,j,k)
      write(ioout, 200) qvec, sqsum(i,j,k)
   enddo
   write(ioout,*)
   enddo
   write(ioout,*)
   enddo
   close(ioout)
   where (nsq.gt.0) sq = sq/dble(nsq)
   !
   outfile = 'sq.dat'
   write(*,'(10x,"Please input the output file name for processed data [sq.dat]: ",$)')
   read(*,'(A)',iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') outfile = trim(oneline)
   !
   open(ioout, file=outfile, status='unknown', action='write')
   write(ioout, 300)
   do i = 2, npt
      if ( nsq(i).lt.1 ) cycle
      qr   = dble(i-1)*dq
      write(ioout, 400) qr, sq(i)
   enddo
   close(ioout)
   !
100 format("# qx",16x,"qy",18x,"qz",18x,"S(q)")
200 format(4G20.10)
300 format("# q          S(q)")
400 format(2G20.10)
return
end subroutine

subroutine outputgr
use sfvars
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i
   real(q) :: r, r2, Rs(2), nnei
   equivalence ( Rs(1), r2 ), ( Rs(2), r )
   !----------------------------------------------------------------------------
   outfile = 'gr.dat'
   write(*,'(/,10x,"Please input the output file name for the calculated g(r)  [gr.dat]: ",$)')
   read(*,'(A)',iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') outfile = trim(oneline)
   !
   open(ioout, file=outfile, status='unknown', action='write')
   write(ioout, 100)
   r    = halfdr
   nnei = 0.D0
   do i = 1, ngr
      r2   = r * r
      nnei = nnei + nnsum(i)
      write(ioout, 200) r, grsum(i)/Rs(ir2r), nnei
      r  = r + dr
   enddo
   close(ioout)
   !
100 format("# r",17x,"g(r)",16x,"Num_of_neighbors" )
200 format(3G20.10)
return
end subroutine

subroutine outputsf_evo
use sfvars
use iounits
implicit none
   !----------------------------------------------------------------------------
   integer :: i
   !----------------------------------------------------------------------------
   outfile = 'sq_time.dat'
   write(*,'(/,10x,"Please input the output file name for the evolution of S(q)[",A,"]: ",$)') trim(outfile)
   read(*,'(A)',iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') outfile = trim(oneline)

   open(ioout, file=outfile, status='unknown', action='write')
   write(ioout, 100) qvec
   do i = istr, iend
      write(ioout, 200) timestep(i), sqsum(i,1,1)
   enddo
   close(ioout)
   !
100 format("# Evolution of S(q) at q=[",3F8.3,"]",/,"# istep   S(q)")
200 format(I10,1x,G20.10)
return
end subroutine
