! Subroutine to analyse the bond type based on HA model
! Only bonded pairs are considered
module HABond
use sfvars
implicit none
   private
   public :: HABond_Init, HABond_Analyse, HABond_Output, HABond_Map
   !----------------------------------------------------------------------------
   integer, allocatable :: int_1d(:), int_2d(:,:)
   integer, allocatable :: nnei(:), nblist(:,:), bonds(:,:)
   integer, allocatable :: bondlist(:,:), typelist(:)
   integer, allocatable :: comnei(:), hasbond(:), profile(:)
   !----------------------------------------------------------------------------
   integer :: i, j, m, n, it, nbondtype, nbonds, maxnei, maxbondtype, bondtype
   integer :: maxbonds
   real(q) :: HArcut, rcutsq
   !----------------------------------------------------------------------------
contains
   !----------------------------------------------------------------------------
   subroutine HABond_Init(flag)
   implicit none
      integer, intent(out):: flag
      !-------------------------------------------------------------------------
      flag = 0
      write(*,'(/,10x,"Please input the cutoff distance for bonds: ", $)')
      read(*,*, iostat=ioerr) HArcut
      if (ioerr.ne.0) then
         flag = 1
         return
      endif
      rcutsq = HArcut * HArcut
      !
      if (allocated(nnei   )) deallocate(nnei  )
      if (allocated(nblist )) deallocate(nblist)
      if (allocated(profile)) deallocate(profile)
      if (allocated(typelist)) deallocate(typelist)
      if (allocated(bonds))    deallocate(bonds)
      maxnei = 20
      maxbondtype = 100
      maxbonds    = nsel*6
      allocate( nnei(nsel), nblist(maxnei, nsel) )
      allocate( profile(maxbondtype), typelist(maxbondtype), bonds(2,maxbonds) )
      profile   = 0
      typelist  = 0
      nbondtype = 0
      !-------------------------------------------------------------------------
   return
   end subroutine HABond_Init
   !
   subroutine HABond_Analyse()
   implicit none
      !-------------------------------------------------------------------------
      integer :: Iimg, II, JJ, KK, LL, flag, bp
      real(q) :: Rij(3), r2, posi(3)
      character(len=1) :: signals(3) = (/ '\', '|', '/' /)
      !-------------------------------------------------------------------------
      write(*,'(10x,"Analysing ...", $)')
      call CPU_TIME( tmbeg )
      do Iimg = istr, iend, inc
         write(*,777) char(8), signals(mod(Iimg,3)+1)
         nnei   = 0
         nblist = 0
         nbonds = 0
         OneImg = atpos(:,:,Iimg)
         prd    = box(:,Iimg)
         !----------------------------------------------------------------------
         do II = 1, nsel-1
            i  = list(II)
            posi = OneImg(:,i)
            !
            do JJ = II+1, nsel
               j   = list(JJ)
               Rij = OneImg(:,j) - posi
               Rij = (Rij - NINT(Rij))
               Rij(1) = Rij(1) * prd(1) + Rij(2)*prd(4) + Rij(3)*prd(5)
               Rij(2) = Rij(2) * prd(2) + Rij(3)*prd(6)
               Rij(3) = Rij(3) * prd(3)
               r2  = sum(Rij*Rij)
               !
               if (r2.le.rcutsq) then
                  nnei(II) = nnei(II) + 1
                  nnei(JJ) = nnei(JJ) + 1
                  nbonds   = nbonds + 1
                  if (nnei(II).gt.maxnei.or.nnei(JJ).gt.maxnei) then
                     if (allocated(int_2d)) deallocate(int_2d)
                     allocate(int_2d(maxnei, nsel))
                     int_2d = nblist
                     deallocate(nblist)
                     maxnei = maxnei + 5
                     allocate(nblist(maxnei, nsel))
                     nblist(1:maxnei-5,:) = int_2d
                  endif
                  nblist(nnei(II),II) = JJ
                  nblist(nnei(JJ),JJ) = II
                  if (nbonds.gt.maxbonds) then
                     if (allocated(int_2d)) deallocate(int_2d)
                     allocate(int_2d(2,maxbonds))
                     int_2d = bonds
                     maxbonds = maxbonds + nsel*2
                     deallocate(bonds)
                     allocate(bonds(2,maxbonds))
                     bonds(:,1:maxbonds-nsel-nsel) = int_2d
                  endif
                  bonds(:,nbonds) = (/ II, JJ /)
               endif
            enddo
         enddo
         !----------------------------------------------------------------------
         ! Now to analyse the bond type
         if (nbonds.lt.1) cycle
         !
         if (allocated(comnei )) deallocate(comnei)
         if (allocated(hasbond)) deallocate(hasbond)
         allocate( comnei(maxnei), hasbond(maxnei) )
         do i = 1, nbonds
            j = 0 ! # of common neighbors
            m = 0 ! # of bonds between common neighbors
            n = 1 ! # extra factor; if all bonds between common neighbors are linked, 2
            comnei  = 0
            hasbond = 0
            II = bonds(1,i)
            JJ = bonds(2,i)
            do KK = 1, nnei(II)
            do LL = 1, nnei(JJ)
               if (nblist(KK, II).eq.nblist(LL, JJ)) then
                  j = j + 1
                  comnei(j) = nblist(KK, II)
               endif
            enddo
            enddo
            do KK = 1, j-1
            do LL = KK+1, j
               flag = bonded(comnei(KK), comnei(LL))
               m = m + flag
               hasbond(KK) = max(hasbond(KK), flag)
               hasbond(LL) = max(hasbond(LL), flag)
            enddo
            enddo
            if (j.eq.4.and.m.eq.2.and.sum(hasbond).lt.4) n = 2
            bondtype = 1000 + j * 100 + m * 10 + n
            bp = lookup(bondtype)
            profile(bp) = profile(bp) + 1
         enddo
         !
      enddo
      call CPU_TIME( tmend )
      write(*, 888) char(8), tmend-tmbeg
777   format(2A1,$)
888 format(A1,"Done!",/,10x,"Total CPU time used: ",F8.2," s")
      !-------------------------------------------------------------------------
   return
   end subroutine HABond_Analyse
   
   integer function bonded(ii, jj)
   implicit none
      !-------------------------------------------------------------------------
      integer, intent(in) :: ii, jj
      integer :: i, j
      !-------------------------------------------------------------------------
      bonded = 0
      do i = 1, nnei(ii)
         if ( jj.eq.nblist(i, ii) ) then
            bonded = 1
            exit
         endif
      enddo
   return
   end function bonded
 
   integer function lookup(ip)
   implicit none
      !-------------------------------------------------------------------------
      integer, intent(in) :: ip
      integer :: i
      !-------------------------------------------------------------------------
      lookup = 0
      do i = 1, nbondtype
         if ( ip.eq.typelist(i) ) lookup = i
      enddo
      if ( lookup.eq.0 ) then
         nbondtype = nbondtype + 1
         if (nbondtype.gt.maxbondtype) then
            if (allocated(int_1d )) deallocate(int_1d)
            allocate(int_1d(maxbondtype))
            int_1d  = typelist
            maxbondtype = maxbondtype + 10
            deallocate( typelist )
            allocate( typelist(maxbondtype) )
            typelist(1:maxbondtype-10) = int_1d
            !
            int_1d = profile
            deallocate( profile )
            allocate( profile(maxbondtype))
            profile(1:maxbondtype-10)  = int_1d
            profile(maxbondtype-10:maxbondtype)  = 0
         endif
         typelist(nbondtype) = ip
         lookup = nbondtype
      endif
      return
      !-------------------------------------------------------------------------
   end function lookup
   !----------------------------------------------------------------------------
   subroutine HABond_Output()
   implicit none
      !-------------------------------------------------------------------------
      real(q)             :: inv_total
      integer             :: i, ntotal, j, idum
      !-------------------------------------------------------------------------
      outfile = 'ha.dat'
      write(*, '(/,10x,"Please input the filename to output profile info[",A,"]: ", $)') trim(outfile)
      read(*, '(A)', iostat=ioerr) oneline
      if (ioerr.eq.0.and.oneline.ne.'') outfile = trim(oneline)
      !
      open(77, file=outfile, action='write', status='unknown')
      !
      ntotal = sum(profile(1:nbondtype))
      inv_total = 1.D0/dble(ntotal)
      !-------------------------------------------------------------------------
      ! sort the bond profile according to # of bonds
      do i = 1, nbondtype-1
      do j = i, nbondtype 
         if ( profile(j).gt.profile(i) ) then
            idum       = profile(j)
            profile(j) = profile(i)
            profile(i) = idum
            idum        = typelist(j)
            typelist(j) = typelist(i)
            typelist(i) = idum
         endif
      enddo
      enddo
      !-------------------------------------------------------------------------
      write(77, '("# HA analyse for frames ", I6," to ", I6, " by ", I6)') istr, iend, inc
      write(77, '("# Average bonds per atom:", F8.3)') dble(ntotal)/dble(nsel*((iend-istr)/inc+1))
      write(77, '("# index  bondtype  num_bonds  ratio")')
      do i = 1, nbondtype
         write(77, 100) i, typelist(i), profile(i), dble(profile(i))*inv_total
      enddo
      close(77)
100   format(I4,1x,I4,1x,I10,2x,F15.6)
      !-------------------------------------------------------------------------
   return
   end subroutine HABond_Output

   subroutine HABond_Map()
   implicit none
      !-------------------------------------------------------------------------
      integer :: Iimg, II, JJ, KK, LL, flag, bp
      real(q) :: Rij(3), r2, posi(3)
      character(len=1) :: signals(3) = (/ '\', '|', '/' /)
      !-------------------------------------------------------------------------
      outfile = 'ha_map.dat'
      write(*, '(/,10x,"Please input the filename to output profile info[",A,"]: ", $)') trim(outfile)
      read(*, '(A)', iostat=ioerr) oneline
      if (ioerr.eq.0.and.oneline.ne.'') outfile = trim(oneline)
      !
      nnei   = 0
      nblist = 0
      nbonds = 0
      OneImg = atpos(:,:,istr)
      prd    = box(:,istr)
      !----------------------------------------------------------------------
      do II = 1, nsel-1
         i  = list(II)
         posi = OneImg(:,i)
         !
         do JJ = II+1, nsel
            j   = list(JJ)
            Rij = OneImg(:,j) - posi
            Rij = (Rij - NINT(Rij))
            Rij(1) = Rij(1) * prd(1) + Rij(2)*prd(4) + Rij(3)*prd(5)
            Rij(2) = Rij(2) * prd(2) + Rij(3)*prd(6)
            Rij(3) = Rij(3) * prd(3)
            r2  = sum(Rij*Rij)
            !
            if (r2.le.rcutsq) then
               nnei(II) = nnei(II) + 1
               nnei(JJ) = nnei(JJ) + 1
               nbonds   = nbonds + 1
               if (nnei(II).gt.maxnei.or.nnei(JJ).gt.maxnei) then
                  if (allocated(int_2d)) deallocate(int_2d)
                  allocate(int_2d(maxnei, nsel))
                  int_2d = nblist
                  deallocate(nblist)
                  maxnei = maxnei + 5
                  allocate(nblist(maxnei, nsel))
                  nblist(1:maxnei-5,:) = int_2d
               endif
               nblist(nnei(II),II) = j
               nblist(nnei(JJ),JJ) = i
               if (nbonds.gt.maxbonds) then
                  if (allocated(int_2d)) deallocate(int_2d)
                  allocate(int_2d(2,maxbonds))
                  int_2d = bonds
                  maxbonds = maxbonds + nsel*2
                  deallocate(bonds)
                  allocate(bonds(2,maxbonds))
                  bonds(:,1:maxbonds-nsel-nsel) = int_2d
               endif
               bonds(:,nbonds) = (/ i, j /)
            endif
         enddo
      enddo
      !----------------------------------------------------------------------
      ! Now to analyse the bond type
      if (nbonds.lt.1) return
      open(77, file=outfile, action='write', status='unknown')
      write(77, '("# Average bonds per atom:", F7.3)') dble(nbonds)/dble(nsel)
      write(77, 199)
      !
      if (allocated(comnei )) deallocate(comnei)
      if (allocated(hasbond)) deallocate(hasbond)
      allocate( comnei(maxnei), hasbond(maxnei) )
      do i = 1, nbonds
         j = 0 ! # of common neighbors
         m = 0 ! # of bonds between common neighbors
         n = 1 ! # extra factor; if all bonds between common neighbors are linked, 2
         comnei  = 0
         hasbond = 0
         II = bonds(1,i)
         JJ = bonds(2,i)
         do KK = 1, nnei(II)
         do LL = 1, nnei(JJ)
            if (nblist(KK, II).eq.nblist(LL, JJ)) then
               j = j + 1
               comnei(j) = nblist(KK, II)
            endif
         enddo
         enddo
         do KK = 1, j-1
         do LL = KK+1, j
            flag = bonded(comnei(KK), comnei(LL))
            m = m + flag
            hasbond(KK) = max(hasbond(KK), flag)
            hasbond(LL) = max(hasbond(LL), flag)
         enddo
         enddo
         if (j.eq.4.and.m.eq.2.and.sum(hasbond).lt.4) n = 2
         bondtype = 1000 + j * 100 + m * 10 + n
         bp = lookup(bondtype)
         !
         Rij = OneImg(:,JJ) - OneImg(:,II)
         Rij = Rij - NINT(Rij)
         posi = ( OneImg(:,II) + 0.5*Rij )
         posi(1) = posi(1) * prd(1) + posi(2)*prd(4) + posi(3)*prd(5)
         posi(2) = posi(2) * prd(2) + posi(3)*prd(6)
         posi(3) = posi(3) * prd(3)
         write(77, 200) i, posi, bp, bondtype
      enddo
      close(77)
      !
199   format("#index x y z bond-type")
200   format(I8,2x,3(F15.8,2x),I3,2x,I4)
      !-------------------------------------------------------------------------
   return
   end subroutine HABond_Map
end module
