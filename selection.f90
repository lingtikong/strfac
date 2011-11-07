! To define the selection of atoms, only info from the first image is used
! selection grammar
!===============================================================================
! Command                              :     selection
!-------------------------------------------------------------------------------
!    all                               :  all atoms
!    type  =  t1 t2 t3 ...             :  select atoms with types t1, t2, t3 and et al.
!    type  != t1 t2 t3 ...             :  not select atoms with types t1, t2, t3 and et al.
!    index >  num                      :  select atoms with index  >  num
!    index >= num                      :  select atoms with index  >= num
!    index <  num                      :  select atoms with index  <  num
!    index =  num1 num2 num3 ...       :  select atoms with index num1, num2, num3 ...
!    index <= num                      :  select atoms with index  <= num
!    index <> num1 num2                :  select atoms with index   num1 <= index <= num2
!    index m= num1 num2                :  select atoms with index  mod(index,num1)= num2
!    x     >= pos                      :  select atoms with fractional coordinates >= pos
!    x     <= pos                      :  select atoms with fractional coordinates <= pos
!    x     <> pos1 pos2                :  select atoms with fractional coordinates within [pos1 pos2]
!-------------------------------------------------------------------------------
!    & and | can be used to combine two or more commands; but the total length of
!    the command cannot be longer than 1024 characters.
!    All position related operations are based on the last image read.
!
! The default of selection is all.
!===============================================================================
subroutine selection
use sfvars
implicit none
   !----------------------------------------------------------------------------
   integer             :: nums(256), ncmd, icmd, i, j, idum, ip, dir
   real(q)             :: fpos(5)
   logical             :: first, inside
   character (len=1024):: command
   character (len=6  ) :: seltype
   character (len=2  ) :: op, logic
   !
   character (len=20), allocatable :: cmdlist(:)
   integer, save       :: nvisit = 0
   !----------------------------------------------------------------------------
   !
   atsel = 1
   !
   do while (.true.)
      write(*, '(/,10x,"Total number of atoms in system is: ", I10)') natom
      write(*, '( 10x,"Please input your selection command, for grammar, type help")')
      write(*, '( 10x,":", $)')
      read(*, '(A)', iostat=ioerr ) command
      if ( ioerr.eq.0.and.(command.eq.'h'.or.command.eq.'help') ) then
         write(*,'(10x,"The available selection commands are:")')
         write(*,'(12x,"all")')
         write(*,'(12x,"type   =  type1 type2 ...")')
         write(*,'(12x,"type   != type1 type2 ...")')
         write(*,'(12x,"index  >  atom_index" )')
         write(*,'(12x,"index  >= atom_index" )')
         write(*,'(12x,"index  <  atom_index" )')
         write(*,'(12x,"index  <= atom_index" )')
         write(*,'(12x,"index  <> atom_index1 atom_index2" )')
         write(*,'(12x,"index  =  num1 num2 ...")')
         write(*,'(12x,"index  m= num1 num2" )')
         write(*,'(12x,"index  m= num0 num1 num2")')
         write(*,'(12x,"x      >  fractiona-position")')
         write(*,'(12x,"x      <  fractiona-position")')
         write(*,'(12x,"x      <> fractiona-position1 fractiona-position2")')
         write(*,'(10x,"You can use & and/or | to make combination.")')
         cycle
      endif
      exit
   enddo
   !
   if ( ioerr.eq.0.and.command.ne.'' ) then
      !
      call parsing  ! to parse the command line
      logic  = '&'
      icmd   = 1
      !
      do while ( icmd.lt.ncmd )
         !
         if ( logic.ne.'&'.and.logic.ne.'|'.and.logic.ne.'&&'.and.logic.ne.'||' ) then
            write(*,'(/,10x,"Error: Unknown logic operation of ", A, "!")') logic
            logic = cmdlist(icmd+1)
            icmd = icmd + 2
            cycle
         endif
         !
         seltype = cmdlist(icmd)
         select case ( seltype )
         case ( 'all' )  ! select all atoms
            select case ( logic )
            case ( '|', '||' )
               atsel = 1
            end select
         case ( 'type' ) ! type operation
            if ( (ncmd-icmd).lt.2 ) exit
            icmd = icmd + 1
            op   = trim(cmdlist(icmd))
            idum = 0
            do i = icmd+1, ncmd
               read(cmdlist(i), *, iostat=ioerr) nums(idum+1)
               if ( ioerr.ne.0 ) exit
               idum = idum + 1
            enddo
            icmd = icmd + idum
            select case ( logic )
            case ( '&', '&&' )
               select case ( op )
               case ( '=' )
                  do i = 1, natom
                     if ( atsel(i).eq.1 ) then
                        ip = attyp(i)
                        inside = .false.
                        do j = 1, idum
                           if (ip.eq.nums(j) ) then
                              inside=.true.
                              exit
                           endif
                        enddo
                        if ( .not.inside ) atsel(i) = 0
                     endif
                  enddo
               case ( '!=' )
                  do i = 1, natom
                     if ( atsel(i).eq.1 ) then
                        ip = attyp(i)
                        inside = .false.
                        do j = 1, idum
                           if (ip.eq.nums(j) ) then
                              inside=.true.
                              exit
                           endif
                        enddo
                        if ( inside ) atsel(i) = 0
                     endif
                  enddo
               case default
                  write(*,'(/,10x,"Error: Unknown type operation of ", A, "!")') op
               end select
            case ( '|', '||' )
               select case ( op )
               case ( '=' )
                  do i = 1, natom
                     if ( atsel(i).ne.1 ) then
                        ip = attyp(i)
                        inside = .false.
                        do j = 1, idum
                           if (ip.eq.nums(j) ) then
                              inside=.true.
                              exit
                           endif
                        enddo
                        if ( inside ) atsel(i) = 1
                     endif
                  enddo
               case ( '!=' )
                 do i = 1, natom
                     if ( atsel(i).ne.1 ) then
                        ip = attyp(i)
                        inside = .false.
                        do j = 1, idum
                           if (ip.eq.nums(j) ) then
                              inside=.true.
                              exit
                           endif
                        enddo
                        if ( .not.inside ) atsel(i) = 1
                     endif
                  enddo
               case default
                  write(*,'(/,10x,"Error: Unknown type operation of ", A, "!")') op
               end select
            end select
         case ( 'index' ) ! index operation
            if ( (ncmd-icmd).lt.2 ) exit
            icmd = icmd + 1
            op = trim(cmdlist(icmd))
            idum = 0
            do i = icmd+1, ncmd
               read(cmdlist(i), *, iostat=ioerr) nums(idum+1)
               if ( ioerr.ne.0 ) exit
               idum = idum + 1
            enddo
            icmd = icmd + idum
            !
            select case (logic)
            case ( '&','&&' )
               select case (op)
               case ( '>'  )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index > operation!")')
                     exit
                  endif
                  if ( nums(1).lt.0.or.nums(1).gt.natom ) then
                     write(*,'(/,10x,"Error: wrong parameter for index > operation!")')
                     exit
                  endif
                  atsel(1:nums(1)) = 0
               case ( '>=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index >= operation!")')
                     exit
                  endif
                  if ( nums(1).lt.1.or.nums(1).gt.natom ) then
                     write(*,'(/,10x,"Error: wrong parameter for index >= operation!")')
                     exit
                  endif
                  atsel(1:nums(1)-1) = 0
               case ( '<'  )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index < operation!")')
                     exit
                  endif
                  if ( nums(1).lt.1.or.nums(1).gt.natom+1 ) then
                     write(*,'(/,10x,"Error: wrong parameter for index < operation!")')
                     exit
                  endif
                  atsel(nums(1):natom) = 0
               case ( '<=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index <= operation!")')
                     exit
                  endif
                  if ( nums(1).lt.1.or.nums(1).gt.natom ) then
                     write(*,'(/,10x,"Error: wrong parameter for index <= operation!")')
                     exit
                  endif
                  atsel(nums(1)+1:natom) = 0
               case ( '='  )
                  do i = 1, natom
                     if ( atsel(i).eq.1 ) then
                        inside = .false.
                        do j = 1, idum
                           if ( i.eq.nums(j) ) then
                              inside = .true.
                              exit
                           endif
                        enddo
                        if ( .not.inside ) atsel(i) = 0
                     endif
                  enddo
               case ( '<>' )
                  if ( idum.ne.2 ) then
                     write(*,'(/,10x,"Error: index <> operation takes and only takes 2 parameters!")')
                     exit
                  endif
                  nums(1) = min(natom, max(1,nums(1)))
                  nums(2) = min(natom, max(1,nums(2)))
                  if ( nums(1).gt.nums(2) ) then
                     write(*,'(/,10x,"Error: wrong arguments for index <> operation!")')
                     exit
                  endif
                  atsel(1:nums(1)-1)     = 0
                  atsel(nums(2)+1:natom) = 0
               case ( 'm=' )
                  if ( idum.lt.2.or.idum.gt.3 ) then
                     write(*,'(/,10x,"Error: index m= operation takes 2-3 parameters!")')
                     exit
                  elseif (idum.eq.3 ) then
                     nums(4) = nums(1)
                     nums(1) = nums(2)
                     nums(2) = nums(3)
                  else
                     nums(4) = 0
                  endif
                  if ( nums(1).lt.1.or.nums(1).gt.natom.or.nums(2).ge.nums(1) ) then
                     write(*,'(/,10x,"Error: wrong arguments for index m= operation!")')
                     exit
                  endif
                  do i = 1, natom
                     if ( atsel(i).eq.1 ) then
                        if ( mod(i-nums(4), nums(1)).ne.nums(2) ) atsel(i) = 0
                     endif
                  enddo
               case default
                  write(*,'(/,10x,"Error: Unknown index operation of ", A, "!")') op
               end select
            case ( '|','||' )
               select case (op)
               case ( '>'  )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index > operation!")')
                     exit
                  endif
                  atsel(nums(1)+1:natom) = 1
               case ( '>=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index >= operation!")')
                     exit
                  endif
                  atsel(nums(1):natom) = 1
               case ( '<'  )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index < operation!")')
                     exit
                  endif
                  atsel(1:nums(1)-1) = 1
               case ( '<=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for index <= operation!")')
                     exit
                  endif
                  atsel(1:nums(1)) = 1
               case ( '='  )
                  do j = 1, idum
                     i = nums(j)
                     if ( i.ge.1.and.i.le.natom) atsel(i) = 1
                  enddo
               case ( '<>' )
                  if ( idum.ne.2 ) then
                     write(*,'(/,10x,"Error: index <> operation takes and only takes 2 parameters!")')
                     exit
                  endif
                  nums(1) = min(natom, max(1,nums(1)))
                  nums(2) = min(natom, max(1,nums(2)))
                  if ( nums(1).gt.nums(2) ) then
                     write(*,'(/,10x,"Error: wrong arguments for index <> operation!")')
                     exit
                  endif
                  atsel(nums(1):nums(2)) = 1
               case ( 'm=' )
                  if ( idum.ne.2.and.idum.ne.3 ) then
                     write(*,'(/,10x,"Error: index m= operation takes 2-3 parameters!")')
                     exit
                  elseif ( idum.eq.3 ) then
                     nums(4) = nums(1)
                     nums(1) = nums(2)
                     nums(2) = nums(3)
                  else
                     nums(4) = 0
                  endif
                  if ( nums(1).lt.1.or.nums(1).gt.natom.or.nums(2).ge.nums(1) ) then
                     write(*,'(/,10x,"Error: wrong arguments for index m= operation!")')
                     exit
                  endif
                  do i = 1, natom
                     if ( mod(i-nums(4), nums(1)).eq.nums(2) ) atsel(i) = 1
                  enddo
               case default
                  write(*,'(/,10x,"Error: Unknown index operation of ", A, "!")') op
               end select
            end select
         case ( 'x', 'y', 'z' ) ! position operation
            if ( (ncmd-icmd).lt.2 ) exit
            icmd = icmd + 1
            op   = trim(cmdlist(icmd))
            idum = 0
            do i = icmd+1, ncmd
               read(cmdlist(i), *, iostat=ioerr) fpos(idum+1)
               if ( ioerr.ne.0 ) exit
               idum = idum + 1
            enddo
            icmd = icmd + idum
            !
            select case ( seltype )
            case ( 'x' )
               dir = 1
            case ( 'y' )
               dir = 2
            case ( 'z' )
               dir = 3
            end select
            !
            select case (logic)
            case ( '&','&&' )
               select case (op)
               case ( '>', '>=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for ",A2," > operation!")') trim(op)
                     exit
                  endif
                  do i = 1, natom
                     if ( OneImg(dir,i).lt.fpos(1) ) atsel(i) = 0
                  enddo
               case ( '<', '<=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for ",A2," < operation!")') trim(op)
                     exit
                  endif
                  do i = 1, natom
                     if ( OneImg(dir,i).gt.fpos(1) ) atsel(i) = 0
                  enddo
               case ( '<>'      )
                  if ( idum.ne.2 ) then
                     write(*,'(/,10x,"Error: ",A2," <> operation takes and only takes 2 parameters!")') trim(op)
                     exit
                  endif
                  do i = 1, natom
                     if ( OneImg(dir,i).lt.fpos(1).or.OneImg(dir,i).gt.fpos(2) ) atsel(i) = 0
                  enddo
               case default
                  write(*,'(/,10x,"Error: Unknown index operation of ", A, "!")') op
               end select
            case ( '|','||' )
               select case (op)
               case ( '>', '>=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for ",A2," > operation!")') trim(op)
                     exit
                  endif
                  do i = 1, natom
                     if ( OneImg(dir,i).ge.fpos(1) ) atsel(i) = 1
                  enddo
               case ( '<', '<=' )
                  if ( idum.gt.1 ) then
                     write(*,'(/,10x,"Error: too many parameters for ",A2," < operation!")') trim(op)
                     exit
                  endif
                  do i = 1, natom
                     if ( OneImg(dir,i).le.fpos(1) ) atsel(i) = 1
                  enddo
               case ( '<>'      )
                  if ( idum.ne.2 ) then
                     write(*,'(/,10x,"Error: ",A2," <> operation takes and only takes 2 parameters!")') trim(op)
                     exit
                  endif
                  do i = 1, natom
                     if ( OneImg(dir,i).ge.fpos(1).and.OneImg(dir,i).le.fpos(2) ) atsel(i) = 1
                  enddo
               case default
                  write(*,'(/,10x,"Error: Unknown index operation of ", A, "!")') op
               end select
            end select
         case default
            write(*,'(10x,"Error: wrong selection keyword of ", A, "!")') trim(seltype)
         end select
         !
         if ( icmd.lt.ncmd ) logic = cmdlist(icmd+1)
         icmd = icmd + 2
      enddo
   endif
   !
   nsel = count(atsel.eq.1)
   write(*,'(10x,"Total number of atoms in your selection is: ", I8)') nsel
   !
   if ( nsel.gt.0 ) then
      if ( allocated(list) ) deallocate( list )
      allocate( list(nsel) )
      list = pack(indice, atsel.eq.1)
   endif
   !----------------------------------------------------------------------------
return
contains
   subroutine parsing
   implicit none
      !-------------------------------------------------------------------------
      integer :: nt, istr, i, j
      !-------------------------------------------------------------------------
      nt = len_trim(command)
      do i = 1, nt
         if ( command(i:i).ne.' ' ) then
            istr = i
            exit
         endif
      enddo
      ncmd = 0
      do i = istr+1, nt
         j = i-1
         if ( command(i:i).eq.' '.and.command(j:j).ne.' ' ) ncmd = ncmd + 1
      enddo
      if ( command(nt:nt).ne.' ' ) ncmd = ncmd + 1
      !
      do i = 1, nt
         if ( command(i:i).ge.'A'.and.command(i:i).le.'Z' ) command(i:i) = char(ichar(command(i:i))+32)
      enddo
      if ( ncmd.ge.1 ) then
         if ( allocated(cmdlist) ) deallocate(cmdlist)
         allocate( cmdlist(ncmd) )
         read(command, *, iostat=ioerr) cmdlist
         if ( ioerr.ne.0 ) then
            write(*,'("Error: error encountered while trying to read selection command!")')
            stop
         endif
      endif
   return
   end subroutine
   !----------------------------------------------------------------------------
end subroutine

subroutine output_selection()
use sfvars
implicit none
   !----------------------------------------------------------------------------
   integer i
   !----------------------------------------------------------------------------
   outfile = 'SelAtom.dat'
   write(*, '(/,10x,"Please input your output file name[",A,"]: ")') trim(outfile)
   read(*,'(A)', iostat=ioerr) oneline
   if (ioerr.eq.0.and.oneline.ne.'') outfile = trim(oneline)
   open(12, file=outfile, action='write', status='unknown')
   write(12,'("# atoms in selection")')
   do i = 1, nsel
      write(12, '(I10)') list(i)
   enddo
   close(12)
   !----------------------------------------------------------------------------
end subroutine
