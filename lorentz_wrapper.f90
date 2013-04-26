! Wrapper to assign C callable names to Fortran routines to be used from 
! Python.
module lorentz_wrapper
use lorentz, only : integrate
use iso_c_binding, only: c_double

implicit none

contains

subroutine c_integrate(data) bind (c)
real(c_double), intent(out) :: data(3, 1000)
call integrate(data)
end subroutine

end module
