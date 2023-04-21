module fortran_struct
  implicit none

  type :: mydataptr
    integer :: counters(5) = 0
  end type mydataptr

contains

  subroutine inc_counter(pdata)
    type(mydataptr), pointer, intent(inout) :: pdata
    integer :: i

    print *, "Received: pdata%counters(1) = ", pdata%counters(1)

    do i = 1, 5
      pdata%counters(i) = 5
    end do
  end subroutine inc_counter
end module fortran_struct