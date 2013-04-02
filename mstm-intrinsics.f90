      module intrinsics
!
! compiler-dependent intrinsic functions.
!
!
!  last revised: 15 January 2011
!

      implicit none
      contains
!
!   system clock
!
         real function mytime()
         implicit none
         real :: etime
         real(4), parameter :: x=0.
         real :: t(2)
!         mytime=secnds(x)
         mytime=etime(t)
         end function mytime

!   flush is one of the best named functions in fortran, although it does not do what I think it should.   This function flushes
!   the buffer to print unit i, so when the unit is flushed, all data written to the open file will appear in the file.
!   Not all compilers have this function.

         subroutine flush(i)
         implicit none
         integer :: i
         flush(i)
         end subroutine flush
!
!   number of command-line arguments.
!
         integer function mstm_nargs()
         implicit none
         integer nargs
!         mstm_nargs=nargs()
         mstm_nargs=iargc()
         end function mstm_nargs
!
!   command line argument retrieval
!
         subroutine mstm_getarg(char)
         implicit none
         integer :: istat
         character(*) :: char
!         call getarg(1,char,istat)
         call getarg(1,char)
         end subroutine mstm_getarg

      end module intrinsics
