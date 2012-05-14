! Subroutine to analyse the Warren-Cowley parameter
! Only bonded pairs are considered
! csro_A = 1 - z_AB / (c_B*Z_A)
module WC_CSRO
use sfvars
implicit none
   private
   public :: CSRO_Init, CSRO_Analyse, CSRO_Output
   !----------------------------------------------------------------------------
   real(q), allocatable :: NumNei(:,:), NatomPerType(:)
   real(q), allocatable :: concetration(:)
   !----------------------------------------------------------------------------
   real(q) :: HArcut, rcutsq
   !----------------------------------------------------------------------------
contains
   !----------------------------------------------------------------------------
   subroutine CSRO_Init()
   implicit none
      !-------------------------------------------------------------------------
      write(*,'(/,10x,"Please input the cutoff distance for bonds: ", $)')
      read(*,*, iostat=ioerr) HArcut
      if (ioerr.ne.0) HArcut = 3.D0
      rcutsq = HArcut * HArcut
      !
      if (allocated(NumNei       )) deallocate(NumNei )
      if (allocated(NatomPerType )) deallocate(NatomPerType )
      if (allocated(concetration )) deallocate(concetration )
      allocate( NumNei(ntype, ntype), NatomPerType(ntype), concetration(ntype) )
      NumNei       = 0.D0
      NatomPerType = 0.D0
      !-------------------------------------------------------------------------
   return
   end subroutine CSRO_Init
   !
   subroutine CSRO_Analyse()
   implicit none
      !-------------------------------------------------------------------------
      integer :: Iimg, II, JJ, i, j, ip, jp, nproc
      real(q) :: Rij(3), r2, posi(3)
      character(len=1) :: signals(4) = (/ '\', '|', '/','-' /)
      !-------------------------------------------------------------------------
      write(*,'(10x,"Analysing ... ", $)')
      call CPU_TIME( tmbeg )
      do Iimg = istr, iend, inc
         write(*,777) char(8), signals(mod(Iimg,4)+1)
         OneImg = atpos(:,:,Iimg)
         !----------------------------------------------------------------------
         do II = 1, nsel
            i  = list(II)
            ip = attyp(i)
            NatomPerType(ip) = NatomPerType(ip) + 1.D0
            !
            posi = OneImg(:,i)
            !
            do JJ = II+1, nsel
               j   = list(JJ)
               jp  = attyp(j)
               Rij = OneImg(:,j) - posi
               Rij = (Rij - NINT(Rij))
               Rij(1) = Rij(1) * prd(1) + Rij(2)*prd(4) + Rij(3)*prd(5)
               Rij(2) = Rij(2) * prd(2) + Rij(3)*prd(6)
               Rij(3) = Rij(3) * prd(3)
               r2  = sum(Rij*Rij)
               !
               if (r2.le.rcutsq) then
                  NumNei(ip,jp) = NumNei(ip,jp) + 1.D0
                  NumNei(jp,ip) = NumNei(jp,ip) + 1.D0
               endif
            enddo
         enddo
         !----------------------------------------------------------------------
      enddo
      concetration = NatomPerType / dble(nsel*((iend-istr)/inc+1))
      call CPU_TIME( tmend )
      write(*, 888) char(8), tmend-tmbeg
777   format(2A1,$)
888 format(A1,"Done!",/,10x,"Total CPU time used: ",F8.2," s")
      !-------------------------------------------------------------------------
   return
   end subroutine CSRO_Analyse
   !----------------------------------------------------------------------------
   subroutine CSRO_Output()
   implicit none
      !-------------------------------------------------------------------------
      real(q)            :: ntotal
      integer            :: i
      character(len=128) :: fmtstr
      !-------------------------------------------------------------------------
      fmtstr = '(/,10x,10("="),??(12("=")))'
      write(fmtstr(16:17), '(I2.2)') ntype
      write(*, fmtstr)
      fmtstr = '(10x,A10,??(5x,I2,5x))'
      write(fmtstr(10:11), '(I2.2)') ntype
      write(*, fmtstr) "Type", (i, i=1,ntype)
      fmtstr = '(10x,10("-"),??(12("-")))'
      write(fmtstr(14:15), '(I2.2)') ntype
      write(*, fmtstr)
      fmtstr = '(10x,A10,??(1x,F10.6,1x))'
      write(fmtstr(10:11), '(I2.2)') ntype
      write(*, fmtstr) "c (%)", concetration * 100.D0
      fmtstr = '(10x,10("="),??(12("=")))'
      write(fmtstr(14:15), '(I2.2)') ntype
      write(*, fmtstr)
      fmtstr = '(10x,A10,??(5x,I2,5x))'
      write(fmtstr(10:11), '(I2.2)') ntype
      write(*, fmtstr) "CSRO", (i, i=1,ntype)
      fmtstr = '(10x,10("-"),??(12("-")))'
      write(fmtstr(14:15), '(I2.2)') ntype
      write(*, fmtstr)
      fmtstr = '(14x,I2,4x,??(1x,F10.6,1x))'
      write(fmtstr(12:13), '(I2.2)') ntype
      do i = 1, ntype
         ntotal = sum(NumNei(i,:))
         write(*, fmtstr) i, 1.D0 -NumNei(i,:)/(ntotal*concetration)
      enddo
      fmtstr = '(10x,10("="),??(12("=")),/)'
      write(fmtstr(14:15), '(I2.2)') ntype
      write(*, fmtstr)
      !-------------------------------------------------------------------------
   return
   end subroutine CSRO_Output
end module
