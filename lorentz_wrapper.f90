module lorentz_wrapper
use lorentz, only : butterfly, assimilate_data
use iso_c_binding, only: c_double

implicit none

contains

subroutine c_butterfly() bind (c)
call butterfly()
end subroutine

subroutine c_assimilate_data(data) bind (c)
real(c_double), intent(out) :: data(279, 4)
call assimilate_data(data)
end subroutine

end module
