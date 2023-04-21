module fortran_client
  use iso_fortran_env, only: dp => real64
  implicit none

  abstract interface
    function sum_callback(a, b) result(res)
      import dp
      real(dp), intent(in) :: a, b
      real(dp) :: res
    end function sum_callback
  end interface

contains
  function compute_media(f, a, b) result(res)
    real(dp), intent(in) :: a, b
    procedure(sum_callback) :: f
    real(dp) :: res

    print *, "Received a: ", a
    print *, "Received b: ", b

    ! Call Julia function to compute res = a + b
    res = f(a, b)/2
  end function compute_media
end module fortran_client