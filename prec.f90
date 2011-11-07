!******************************************************************
!***  Module of prec, to adapt double precision to all machine. ***
!******************************************************************
module prec
implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
end module

module iounits
implicit none
   integer, parameter :: ioin  = 8
   integer, parameter :: ioout = 9
   integer, parameter :: iotmp = 10
end module
