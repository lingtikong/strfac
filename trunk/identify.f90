! Subroutine to identify the total number of atoms, atom types and first box size
subroutine identify
use sfvars
use iounits
use errmesg
implicit none
   !----------------------------------------------------------------------------
   logical :: fexst
   integer :: i, id, it, nlines, rderr(5), ntm, istep, inttmp(3)
   real(q) :: dptmp(3)
   integer :: nwarning = 0
   !----------------------------------------------------------------------------
   real(q), allocatable :: postmp(:,:,:), boxtmp(:,:), tstmp(:)
   integer, allocatable :: wraptmp(:,:,:)
   real(q), parameter   :: zero = 1.D-6
   character(len=1)     :: signals(4) = (/ '\', '|', '/','-' /)
   character(len=256)   :: fmtstr
   !----------------------------------------------------------------------------
   subname = 'identify'
   !
   call CPU_TIME( tmbeg )
   !
   ! Check if input file exists
   inquire( file=infile, exist=fexst )
   if ( .not.fexst ) then
      write(errinfo, '("File:",A," not found!")') trim(infile)
      call error(1)
   endif
   write(*,'(/, 10x,"Now to identify the dump file, wait for a while ... ", $)')
   !
   ! Now to read the first image
   open(ioin, file=trim(infile), status='old', action='read', iostat=ioerr )
   read(ioin, '(A)', iostat=rderr(1) ) oneline
   read(ioin, *,     iostat=rderr(2) ) istep
   read(ioin, '(A)', iostat=rderr(3) ) oneline
   if ( any(rderr(1:3).ne.0) ) then
      write(errinfo, '("Error encountered while reading lammps atom dump file: ", A)') trim(infile)
      call error(1)
   endif
   read(ioin, *, iostat=ioerr) natom
   if ( ioerr.ne.0 ) then
      write(errinfo, '("Error while reading # of atoms from lammps atom dump file: ", A)') trim(infile)
      call error(ioerr)
   endif
   xy = 0.D0
   xz = 0.D0
   yz = 0.D0
   read(ioin, '(A)', iostat=rderr(1) ) oneline
   read(ioin, '(A)', iostat=rderr(1) ) oneline
   read(oneline,*,iostat=ioerr) xlo, xhi, xy
   if (ioerr.ne.0) read(oneline,*,iostat=rderr(2)) xlo, xhi

   read(ioin, '(A)', iostat=rderr(1) ) oneline
   read(oneline,*,iostat=ioerr) ylo, yhi, xz
   if (ioerr.ne.0) read(oneline,*,iostat=rderr(3)) ylo, yhi

   read(ioin, '(A)', iostat=rderr(1) ) oneline
   read(oneline,*,iostat=ioerr) zlo, zhi, yz
   if (ioerr.ne.0) read(oneline,*,iostat=rderr(4)) zlo, zhi

   read(ioin, '(A)', iostat=rderr(5) ) oneline
   lwrapinfo = len_trim(oneline).ge.37
   if ( any(rderr(1:5).ne.0) ) then
      write(errinfo, '("Error while reading box info from lammps atom dump file: ", A)') trim(infile)
      call error(2)
   endif
   if ( natom.lt.1 ) then
      write(errinfo, '("No atom was identified from lammps atom dump file: ", A)') trim(infile)
      call error(1)
   endif

   xlo = xlo - MIN(0.D0, MIN(xy, MIN(xz, xy+xz)));
   xhi = xhi - MAX(0.D0, MAX(xy, MAX(xz, xy+xz)));

   ylo = ylo - MIN(0.D0,yz);
   yhi = yhi - MAX(0.D0,yz);

   Lx = xhi - xlo
   Ly = yhi - ylo
   Lz = zhi - zlo
   vol = Lx * Ly * Lz
   !
   nwarning = 0
   if ( xy*xy + xz*xz + yz*yz > 1.e-8 ) nwarning = 1
   !
   NimgMax = 500
   ntype   = 0
   !
   if ( allocated( atpos ) ) deallocate( atpos )
   if ( allocated( attyp ) ) deallocate( attyp )
   if ( allocated( atsel ) ) deallocate( atsel )
   if ( allocated( indice) ) deallocate( indice)
   if ( allocated( OneImg) ) deallocate( OneImg)
   if ( allocated( box   ) ) deallocate( box   )
   if ( allocated(timestep)) deallocate(timestep)
   if ( allocated( wrap   )) deallocate( wrap   )
   allocate( atpos(3,natom,NimgMax), box(6,NimgMax), attyp(natom), atsel(natom), indice(natom), OneImg(3,natom), timestep(NimgMax) )
   if ( lwrapinfo ) allocate( wrap(3,natom, NimgMax), wrapone(3,natom) )
   forall (i=1:natom) indice(i) = i
   !
   if ( lwrapinfo ) then
      do i = 1, natom
         read(ioin, *, iostat=ioerr) id, it, dptmp, inttmp
         if ( ioerr.ne.0 ) then
            write(errinfo, '("Error while identifing atom types from lammps atom dump file: ", A)') trim(infile)
            call error(ioerr)
         endif
         attyp(id)     = it
         OneImg(:,id)  = dptmp
         wrapone(:,id) = inttmp
         call findtype(it)
      enddo
   else
      do i = 1, natom
         read(ioin, *, iostat=ioerr) id, it, dptmp
         if ( ioerr.ne.0 ) then
            write(errinfo, '("Error while identifing atom types from lammps atom dump file: ", A)') trim(infile)
            call error(ioerr)
         endif
         attyp(id)   = it
         OneImg(:,id) = dptmp
         call findtype(it)
      enddo
   endif
   !
   nimage = 1
   atpos(:,:,nimage) = OneImg
   timestep(nimage)  = istep
   box(:,1) = (/ Lx, Ly, Lz, xy, xz, yz /)
   !
   read(ioin, '(A)', iostat=ioerr ) oneline
   do while ( ioerr.eq.0 )
      read(ioin, *,     iostat=rderr(1) ) istep
      read(ioin, '(A)', iostat=rderr(2) ) oneline
      if ( any(rderr(1:2).ne.0) ) exit
      !
      read(ioin, *, iostat=ioerr) ntm
      if ( ioerr.ne.0.or.ntm.ne.natom ) exit
      !
      xy = 0.D0
      xz = 0.D0
      yz = 0.D0
      !
      read(ioin, '(A)', iostat=rderr(1) ) oneline
      read(ioin, '(A)', iostat=rderr(2) ) oneline
      read(oneline,*,iostat=ioerr) xlo, xhi, xy
      if (ioerr.ne.0) read(oneline,*,iostat=rderr(2)) xlo, xhi
   
      read(ioin, '(A)', iostat=rderr(3) ) oneline
      read(oneline,*,iostat=ioerr) ylo, yhi, xz
      if (ioerr.ne.0) read(oneline,*,iostat=rderr(3)) ylo, yhi
   
      read(ioin, '(A)', iostat=rderr(4) ) oneline
      read(oneline,*,iostat=ioerr) zlo, zhi, yz
      if (ioerr.ne.0) read(oneline,*,iostat=rderr(4)) zlo, zhi

      read(ioin, '(A)', iostat=rderr(5) ) oneline
      if ( any(rderr(1:5).ne.0) ) exit

      xlo = xlo - MIN(0.D0, MIN(xy, MIN(xz, xy+xz)));
      xhi = xhi - MAX(0.D0, MAX(xy, MAX(xz, xy+xz)));

      ylo = ylo - MIN(0.D0,yz);
      yhi = yhi - MAX(0.D0,yz);

      if ( xy*xy + xz*xz + yz*yz > 1.e-8 ) nwarning = nwarning + 1
      Lx = xhi - xlo
      Ly = yhi - ylo
      Lz = zhi - zlo
      !
      if ( lwrapinfo ) then
         do i = 1, natom
            read(ioin, *, iostat=ioerr) id, it, dptmp, inttmp
            if ( ioerr.ne.0.or.it.ne.attyp(id) ) exit
            OneImg(:,id)  = dptmp
            wrapone(:,id) = inttmp
         enddo
      else
         do i = 1, natom
            read(ioin, *, iostat=ioerr) id, it, dptmp
            if ( ioerr.ne.0.or.it.ne.attyp(id) ) exit
            OneImg(:,id) = dptmp
         enddo
      endif
      if ( i.lt.natom ) exit
      !
      if ( nimage.eq.NimgMax ) then
         if ( allocated( postmp ) ) deallocate( postmp )
         if ( allocated( boxtmp ) ) deallocate( boxtmp )
         if ( allocated( tstmp  ) ) deallocate( tstmp  )
         if ( allocated( wraptmp) ) deallocate( wraptmp)
         allocate( postmp(3, natom, nimage), boxtmp(6,nimage), tstmp(nimage) )
         if (lwrapinfo) then
            allocate(wraptmp(3,natom, nimage))
            wraptmp = wrap
         endif
         postmp  = atpos
         boxtmp  = box
         tstmp   = timestep
         
         NimgMax = NimgMax + 500
         deallocate( atpos, box, timestep )
         allocate( atpos(3, natom, NimgMax), box(6,NimgMax), timestep(NimgMax) )
         if (lwrapinfo) then
            deallocate(wrap)
            allocate( wrap(3,natom,NimgMax) )
            wrap(:,:,1:nimage) = wraptmp
         endif
         atpos(:,:,1:nimage) = postmp
         box(:,1:nimage) = boxtmp
         timestep(1:nimage) = tstmp
      endif
      nimage = nimage + 1
      atpos(:,:,nimage) = OneImg
      if (lwrapinfo) wrap(:,:,nimage) = wrapone
      box(:,nimage) = (/ Lx, Ly, Lz, xy, xz, yz /)
      timestep(nimage) = istep
      !
      it = mod(nimage,4) + 1
      write(*,'(2A1,$)') char(8), signals(it)
      !
      read(ioin, '(A)', iostat=ioerr ) oneline
   enddo
   close(ioin)
   if ( allocated( postmp ) ) deallocate( postmp )
   if ( allocated( boxtmp ) ) deallocate( boxtmp )
   if ( allocated( tstmp  ) ) deallocate( tstmp  )
   if ( allocated( wraptmp) ) deallocate( wraptmp)
   !
   write(*,'(A1,A)') char(8), 'Done!'
   !
   fmtstr   = '(/,10x,"Identified ",I??," images from file ",A,": ", I??," atoms with ", I??, " types. ")'
   rderr(1) = int(log10(dble(nimage))) + 1
   rderr(2) = int(log10(dble(natom ))) + 1
   rderr(3) = int(log10(dble(ntype ))) + 1
   write(fmtstr(23:24),'(I2.2)') rderr(1)
   write(fmtstr(56:57),'(I2.2)') rderr(2)
   write(fmtstr(76:77),'(I2.2)') rderr(3)
   !
   write(*, fmtstr) nimage, trim(infile), natom, ntype
   fmtstr = '(12x,"There are ", I??, " atoms of type ",I??)'
   write(fmtstr(21:22),'(I2.2)') rderr(2)
   write(fmtstr(44:45),'(I2.2)') rderr(3)
   do i = 1, ntype
      write(*, fmtstr) count(attyp.eq.i), i
   enddo
   if (nwarning > 0) then
      fmtstr = '(/,10x,"WARNNING: ", I??," of ", I??, " images are non-orthogonal!! Some results might be unreliable!")'
      write(fmtstr(23:24),'(I2.2)') int(log10(dble(nwarning)))+1
      write(fmtstr(35:36),'(I2.2)') rderr(1)
      write(*, fmtstr) nwarning, nimage
   endif
   !
   call CPU_TIME( tmend )
   tmread = tmend - tmbeg
   write(*,'(10x,"Total CPU time used to read file:",F7.2," s")') tmread
   !----------------------------------------------------------------------------   
return
contains
   subroutine findtype( new )
   implicit none
      !-------------------------------------------------------------------------
      integer, intent(in) :: new
      integer, save :: founds(max_type)
      integer       :: i, j, ip
      !-------------------------------------------------------------------------
      ip = 0
      do i = 1, ntype
         if ( new.eq.founds(i) ) then
            ip = i
            exit
         endif
      enddo
      if ( ip.eq.0 ) then
         ntype = ntype + 1
         if ( ntype.gt.max_type ) stop 'Too many types found!'
         founds(ntype) = new
      endif
      !-------------------------------------------------------------------------
   return
   end subroutine findtype
end subroutine
