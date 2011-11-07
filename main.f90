! Program to evaluate the structure factor for certain atoms from a lammps atom style dump file
program main
use prec
use sfvars
use iounits
use errmesg
use HABond
use WC_CSRO
use SFDist
implicit none
   !----------------------------------------------------------------------------
   integer          :: idum, narg, newselect
   character(len=5) :: jobstr(12) = (/ "S(q) ","Sq(t)", "Sq(r)", "g(r) ", "g(r) ", &
                     & "g(r) ", "g(r) ", "HA   ", "HA   ", "CSRO ", "dump ", "S(t) " /)
   character(len=1) :: signals(3) = (/ '\', '|', '/' /)
   !----------------------------------------------------------------------------
   ! Get dump file name
   narg = iargc()
   if ( narg.gt.0 ) then
      call getarg(1, infile)
   else
      write(*,'(/,10x,"Please input the name of LAMMPS atom dump file: ",$)')
      read(*, *, iostat=ioerr) infile
   endif
   !
   ! Identify # atom, atom types, box sizes, and total # of images
   call identify
   !
   do while (.true.)
      !
      ! Select atoms to evaluate S(q) and/or g(r)
      OneImg = atpos(:,:,nimage)
      nsel = 0
      do while ( nsel.lt.1 )
         call selection
      enddo
      !
      do while ( .true. )
         ! ask job type
         write(*,'(/,10x,"Please select job type:")')
         write(*,'(10x,"=============================================")')
         write(*,'(12x,"1. Structure factor calculation;")')
         write(*,'(12x,"2. Structure factor evolution;")')
         write(*,'(12x,"3. Structure factor distribution;")')
         write(*,'(10x,"---------------------------------------------")')
         write(*,'(12x,"4. Pair correlation calculation;")')
         write(*,'(12x,"5. Planar pair-correlation calculation;")')
         write(*,'(12x,"6. Pair correlation of desired type pair;")')
         write(*,'(12x,"7. Pair correlation of desired groups;")')
         write(*,'(10x,"---------------------------------------------")')
         write(*,'(12x,"8. HA bond analyse;")')
         write(*,'(12x,"9. HA bond map of one image;")')
         write(*,'(11x,"10. Warren-Cowley CSRO parameter;")')
         write(*,'(10x,"---------------------------------------------")')
         write(*,'(11x,"11. dump each selected atom as a single frame")')
         write(*,'(11x,"12. Self-intermediate scattering function")')
         write(*,'(11x,"13. Output indices of atoms in selection")')
         write(*,'(10x,"---------------------------------------------")')
         write(*,'(11x,"14. New selection;")')
         newselect = 14
         write(*,'(12x,"0. Exit.")')
         write(*,'(10x,"Your choice [0]: ", $)')
         read(*,'(A)',iostat= ioerr) oneline
         write(*,'(10x,"=============================================")')
         job = 0
         if ( ioerr.eq.0.and.oneline.ne.'' ) then
            read(oneline, *, iostat=ioerr) job
            if ( ioerr.ne.0 ) job = 0
         endif
         if ( job.lt.1.or.job.gt.newselect ) then
            stop
         elseif ( job.eq.13) then
            call output_selection
            cycle
         elseif ( job.eq.newselect ) then
            exit
         endif
         !
         ! Get the image range to do the real job
         write(*, 100 ) nimage
         write(*, 200 ) jobstr(job), nimage
100      format(/,10x,"Total # of images in the dump file is ", I7 )
200      format(10x,"Please input the desired image range to evaluate ",A4,"[1:",I7,"]: ",$)
201      format(10x,"Your desired image range is [",I7," -",I7,"] with incremental of", I7,".")
         !
         inc = 1
         read(*, '(A)', iostat=ioerr) oneline
         if ( ioerr.ne.0 ) cycle
         if ( oneline.eq.'' ) then
            istr = 1
            iend = nimage
         else
            read(oneline, *, iostat=ioerr) istr, iend, inc
            if ( ioerr.ne.0 ) then
               read(oneline, *, iostat=ioerr) istr, iend
               if ( ioerr.ne.0 ) cycle
               inc = 1
            endif
            if ( iend.lt.istr.or.inc.gt.iend ) cycle
            if ( istr.lt.1    ) istr = 1
            if ( iend.gt.nimage ) iend = nimage
         endif
         write(*,201) istr, iend, inc
         !
         nimgused = (iend-istr)/inc + 1
         !
         ! Initialization
         select case ( job )
         case ( 1 ) ! To calculate S(q)
            !
            write(*,'(/,10x,"Please input the size of q-mesh: ",$)')
            read(*,*,iostat=ioerr) Nx, Ny, Nz
            if ( ioerr.ne.0.or.Nx.lt.1.or.Ny.lt.1.or.Nz.lt.1 ) then
               write(errinfo, '("Wrong input!")')
               call warn(ioerr)
               cycle
            endif
            !
            call CPU_TIME( tmbeg )
            write(*, 999)
            !
            if ( allocated(sqsum) ) deallocate( sqsum )
            if ( allocated(sqcur) ) deallocate( sqcur )
            if ( allocated(qxrx ) ) deallocate( qxrx  )
            if ( allocated(qyry ) ) deallocate( qyry  )
            if ( allocated(qzrz ) ) deallocate( qzrz  )
            allocate( sqsum(Nx,Ny,Nz), sqcur(Nx,Ny,Nz), qxrx(nsel,Nx), qyry(nsel,Ny), qzrz(nsel,Nz) )
            !
            sqsum = 0.D0
            do image = istr, iend, inc
               OneImg = atpos(:,:,image)
               call strfaccal
               sqsum = sqsum + sqcur *conjg(sqcur)
               write(*, 777) char(8), signals(mod(image,3)+1)
            enddo
            sqsum = sqsum/dble(nimgused)
            
            call CPU_TIME( tmend )
            write(*, 888) char(8), tmend-tmbeg
            !
            call outputsf
            !
         case ( 2 ) ! To calculate the evolution of S(q) at a specific q
            write(*,'(/,10x,"Please input the vector q in units of (1/Lx 1/Ly 1/Lz): ",$)')
            read(*,*,iostat=ioerr) qvec
            if ( ioerr.ne.0 ) then
               write(errinfo, '("Wrong input!")')
               call warn(ioerr)
               cycle
            endif
            !
            call CPU_TIME( tmbeg )
            write(*, 999)
            !
            if ( allocated(sqsum) ) deallocate( sqsum )
            if ( allocated(sqcur) ) deallocate( sqcur )
            allocate( sqsum(istr:iend,1,1), sqcur(istr:iend,1,1) )
            !
            do image = istr, iend, inc
               OneImg = atpos(:,:,image)
               call strfacOne
               sqsum(image,1,1) = sqcur(image,1,1) *conjg(sqcur(image,1,1))
               write(*, 777) char(8), signals(mod(image,3)+1)
            enddo
            sqsum = sqsum/dble(nsel)
            !
            call CPU_TIME( tmend )
            write(*, 888) char(8), tmend-tmbeg
            !
            call outputsf_evo
            !
         case ( 3 ) ! Sq(r)
            call StrFacDist()
         case ( 4 ) ! To calculate g(r)
            !
            ir2r = 1
            rmax = 0.5D0*minval( box(1:3,1:nimage) )
            write(*, 300) rmax
300         format(/,10x,"Please input the max value of r in g(r) [",F6.3,"] : ", $)
            read(*,'(A)',iostat=ioerr) oneline
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               dr = rmax
               read(oneline,*,iostat=ioerr) rmax
               if (ioerr.ne.0.or.rmax.gt.dr) rmax = dr
            endif
            !
            write(*,'(10x,"Please input the total # of g(r) points [100]: ", $)')
            read(*,'(A)',iostat=ioerr) oneline
            ngr = 100
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               read(oneline,*,iostat=ioerr) ngr
               if (ioerr.ne.0) ngr = 100
            endif
            if ( ngr.lt.1 ) then
               write(errinfo, '("Number of g(r) points must be greater than 1!")')
               call warn(1-ngr)
               cycle
            endif
            !
            call CPU_TIME( tmbeg )
            write(*, 999)
            !
            plane  = 1
            dr     = rmax/dble(ngr)
            rdr    = 1.D0/dr
            halfdr = dr * 0.5D0
            grbase = 1.D0/(tpi*dr*dble(nsel)*dble(nsel)*dble(nimgused)) ! volume will be multiplied for each image.
            !
            if ( allocated( grsum ) ) deallocate( grsum )
            if ( allocated( nnsum ) ) deallocate( nnsum )
            if ( allocated( nncur ) ) deallocate( nncur )
            allocate( grsum(ngr), nnsum(ngr), nncur(ngr) )
            !
            grsum = 0.D0
            nnsum = 0.D0
            do image = istr, iend, inc
               prd   = box(:,image)
               vol   = prd(1)*prd(2)*prd(3)
               grwt  = grbase * vol
               OneImg = atpos(:,:,image)
               !
               call paircorr
               !
               grsum = grsum + dble(nncur)*grwt
               nnsum = nnsum + dble(nncur)
               !
               write(*, 777) char(8), signals(mod(image,3)+1)
            enddo
            nnsum = nnsum / ( dble(nsel)*dble(nimgused) )
            !
            call CPU_TIME( tmend )
            write(*, 888) char(8), tmend-tmbeg
            !
            call outputgr
            !
         case ( 5 ) ! To calculate planar g(r)
            ir2r = 2
            !
            rmax = 0.5D0*minval( box(1:3,1:nimage) )
            write(*, 300) rmax
            read(*,'(A)',iostat=ioerr) oneline
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               dr = rmax
               read(oneline,*,iostat=ioerr) rmax
               if (ioerr.ne.0.or.rmax.gt.dr) rmax = dr
            endif
            !
            write(*,'(10x,"Please input the total # of g(r) points [100]: ", $)')
            read(*,'(A)',iostat=ioerr) oneline
            ngr = 100
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               read(oneline,*,iostat=ioerr) ngr
               if (ioerr.ne.0) ngr = 100
            endif
            if ( ngr.lt.1 ) then
               write(errinfo, '("Number of g(r) points must be greater than 1!")')
               call warn(1-ngr)
               cycle
            endif
            !
            plane  = 1
            write(*,'(/,10x,"Please select the plane: ")')
            write(*,'(12x,"1. xy;")')
            write(*,'(12x,"2. xz;")')
            write(*,'(12x,"3. yz;")')
            write(*,'(10x,"Your choice [1]: ",$)')
            read(*,'(A)', iostat=ioerr ) oneline
            idum = 1
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               read(oneline,*,iostat=ioerr) idum
               if ( ioerr.ne.0 ) idum = 1
            endif
            !
            call CPU_TIME( tmbeg )
            write(*,999)
            !
            idum = min(3,max(1,idum))
            plane(4-idum) = 0
            !
            dr     = rmax/dble(ngr)
            rdr    = 1.D0/dr
            halfdr = dr * 0.5D0
            grbase = 1.D0/(tpi*dr*dble(nsel)*dble(nsel)*dble(nimgused)) ! volume will be multiplied for each image.
            !
            if ( allocated( grsum ) ) deallocate( grsum )
            if ( allocated( nnsum ) ) deallocate( nnsum )
            if ( allocated( nncur ) ) deallocate( nncur )
            allocate( grsum(ngr), nnsum(ngr), nncur(ngr) )
            !
            grsum = 0.D0
            nnsum = 0.D0
            do image = istr, iend, inc
               prd   = box(:,image)
               vol   = prd(1)*prd(2)*prd(3)
               grwt  = grbase * vol * 2.D0 / sum(prd(1:3)*dble(1-plane))
               prd(1:3) = prd(1:3) * dble(plane)
               OneImg = atpos(:,:,image)
               !
               call paircorr
               !
               grsum = grsum + dble(nncur)*grwt
               nnsum = nnsum + dble(nncur)
               write(*,777) char(8), signals(mod(image,3)+1)
            enddo
            nnsum = nnsum / (dble(nsel)*dble(nimgused))
            !
            call CPU_TIME( tmend )
            write(*, 888) char(8), tmend- tmbeg
            !
            call outputgr
            !
         case ( 6 ) ! To calculate the pair-correlation of desired pair
            ir2r = 1
            !
            rmax = 0.5D0*minval( box(1:3,1:nimage) )
            write(*, 300) rmax
!400         format(/,10x,"Please input the max value of r in g(r) [",F6.3,"] : ", $)
            read(*,'(A)',iostat=ioerr) oneline
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               dr = rmax
               read(oneline,*,iostat=ioerr)  rmax
               if (ioerr.ne.0.or.rmax.gt.dr) rmax = dr
            endif
            !
            write(*,'(10x,"Please input the total # of g(r) points [100]: ", $)')
            read(*,'(A)',iostat=ioerr) oneline
            ngr = 100
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               read(oneline,*,iostat=ioerr) ngr
               if (ioerr.ne.0) ngr = 100
            endif
            if ( ngr.lt.1 ) then
               write(errinfo, '("Number of g(r) points must be greater than 1!")')
               call warn(1-ngr)
               cycle
            endif
            !
            write(*,'(10x,"Please input the type index of the source,   0 for all: ", $)')
            read(*,*,iostat=ioerr) typsrc
            if (ioerr.ne.0.or.typsrc.lt.0.or.typsrc.gt.ntype) cycle
            write(*,'(10x,"Please input the type index of the neighbor, 0 for all: ", $)')
            read(*,*,iostat=ioerr) typdes
            if (ioerr.ne.0.or.typdes.lt.0.or.typdes.gt.ntype) cycle
            !
            call CPU_TIME( tmbeg )
            write(*, 999)
            !
            plane  = 1
            dr     = rmax/dble(ngr)
            rdr    = 1.D0/dr
            halfdr = dr * 0.5D0
            if (typsrc.eq.0) then
               nsrc = nsel
            else
               nsrc   = count(attyp.eq.typsrc.and.atsel.eq.1)
            endif
            if (typdes.eq.0) then
               ndes = nsel
            else
               ndes   = count(attyp.eq.typdes.and.atsel.eq.1)
            endif
            if (nsrc.lt.1.or.ndes.lt.1) cycle
            !
            if (allocated(srclist)) deallocate(srclist)
            if (allocated(deslist)) deallocate(deslist)
            allocate( srclist(nsrc), deslist(ndes) )
            if (typsrc.eq.0) then
               srclist = pack(indice, atsel.eq.1)
            else
               srclist = pack(indice, attyp.eq.typsrc.and.atsel.eq.1)
            endif
            if (typdes.eq.0) then
               deslist = pack(indice, atsel.eq.1)
            else
               deslist = pack(indice, attyp.eq.typdes.and.atsel.eq.1)
            endif
            grbase = 1.D0/(2.D0*tpi*dr*dble(nsrc)*dble(ndes)*dble(nimgused)) ! volume will be multiplied for each image.
            !
            if ( allocated( grsum ) ) deallocate( grsum )
            if ( allocated( nnsum ) ) deallocate( nnsum )
            if ( allocated( nncur ) ) deallocate( nncur )
            allocate( grsum(ngr), nnsum(ngr), nncur(ngr) )
            !
            grsum = 0.D0
            nnsum = 0.D0
            do image = istr, iend, inc
               prd   = box(:,image)
               vol   = prd(1)*prd(2)*prd(3)
               grwt  = grbase * vol
               OneImg = atpos(:,:,image)
               !
               call pppc
               !
               grsum = grsum + dble(nncur)*grwt
               nnsum = nnsum + dble(nncur)
               !
               write(*, 777) char(8), signals(mod(image,3)+1)
            enddo
            nnsum = nnsum / ( dble(nsrc)*dble(nimgused) )
            !
            call CPU_TIME( tmend )
            write(*, 888) char(8), tmend-tmbeg
            !
            call outputgr
            deallocate(srclist, deslist)
            !
         case ( 7 ) ! To calculate the pair-correlation of desired groups
            ir2r = 1
            !
            rmax = 0.5D0*minval( box(1:3,1:nimage) )
            write(*, 300) rmax
            read(*,'(A)',iostat=ioerr) oneline
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               dr = rmax
               read(oneline,*,iostat=ioerr)  rmax
               if (ioerr.ne.0.or.rmax.gt.dr) rmax = dr
            endif
            !
            write(*,'(10x,"Please input the total # of g(r) points [100]: ", $)')
            read(*,'(A)',iostat=ioerr) oneline
            ngr = 100
            if ( ioerr.eq.0.and.oneline.ne.'' ) then
               read(oneline,*,iostat=ioerr) ngr
               if (ioerr.ne.0) ngr = 100
            endif
            if ( ngr.lt.1 ) then
               write(errinfo, '("Number of g(r) points must be greater than 1!")')
               call warn(1-ngr)
               cycle
            endif
            !
            write(*,'(//,10x,"Please select the atoms server as the source:")')
            OneImg = atpos(:,:,nimage)
            call selection
            if (nsel < 1) cycle
            nsrc = nsel
            if (allocated(srclist)) deallocate(srclist)
            allocate( srclist(nsrc) )
            srclist = list
            
            write(*,'(//,10x,"Please select the atoms server as the neighbors:")')
            call selection
            if (nsel < 1) cycle
            ndes = nsel
            if (allocated(deslist)) deallocate(deslist)
            allocate( deslist(ndes) )
            deslist = list
            !
            call CPU_TIME( tmbeg )
            write(*, 999)
            !
            plane  = 1
            dr     = rmax/dble(ngr)
            rdr    = 1.D0/dr
            halfdr = dr * 0.5D0
            !
            grbase = 1.D0/(2.D0*tpi*dr*dble(nsrc)*dble(ndes)*dble(nimgused)) ! volume will be multiplied for each image.
            !
            if ( allocated( grsum ) ) deallocate( grsum )
            if ( allocated( nnsum ) ) deallocate( nnsum )
            if ( allocated( nncur ) ) deallocate( nncur )
            allocate( grsum(ngr), nnsum(ngr), nncur(ngr) )
            !
            grsum = 0.D0
            nnsum = 0.D0
            do image = istr, iend, inc
               prd   = box(:,image)
               vol   = prd(1)*prd(2)*prd(3)
               grwt  = grbase * vol
               OneImg = atpos(:,:,image)
               !
               call pppc
               !
               grsum = grsum + dble(nncur)*grwt
               nnsum = nnsum + dble(nncur)
               !
               write(*, 777) char(8), signals(mod(image,3)+1)
            enddo
            nnsum = nnsum / ( dble(nsrc)*dble(nimgused) )
            !
            call CPU_TIME( tmend )
            write(*, 888) char(8), tmend-tmbeg
            !
            call outputgr
            deallocate(srclist, deslist)
            !
         case ( 8 ) ! HA bond analyse
            call HABond_Init(idum)
            if (idum.eq.0) then
               call HABond_Analyse()
               call HABond_Output()
            endif
         case ( 9 ) ! HA bond map of one image
            call HABond_Init(idum)
            if (idum.eq.0) call HABond_Map()
         case ( 10 ) ! Warren-Cowley chemical short range order parameter
            call CSRO_Init()
            call CSRO_Analyse()
            call CSRO_Output()
         case ( 11 ) ! Write each selected atom as a seperate frame
            call separate()
         case ( 12 ) ! to calculate the self-intermediate scattering function
            if (lwrapinfo) then
               call SISF()
            else
               write(*,'(10x,"Wrap info not available, computation not possible!")')
            endif
         end select
      enddo
   enddo
   !
999 format(10x,"Computing ... ",$)
777 format(2A1,$)
888 format(A1,"Done!",/,10x,"Total CPU time used: ",F8.2," s")
   !
   call free_all
   !----------------------------------------------------------------------------
end program
