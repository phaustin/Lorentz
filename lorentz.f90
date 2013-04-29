! Module to integrate the Lorentz equations; the integrate() subroutine is
! assigned a C callable name in lorentz_wrapper.f90 for use in Python. 
module lorentz
use types, only: dp
use integrators, only: fwd_euler, midpoint, rk4, heun, mod_euler, ck4, ck5
implicit none

contains

subroutine integrate(ys)
! Integration parameters
 ! Integration parameters
real(dp), parameter :: dt = 0.01_dp
integer, parameter :: steps = 1000
integer :: i
! Integration variables
real(dp) :: y(3)
! Lorentz attractor
real(dp), intent(out) :: ys(3, steps)

! Initial conditions
y = [13.0_dp, 8.1_dp, 45.0_dp]

! Integrator loop
do i = 0, steps-1
    ! Change integrator here as needed
    y = ck4(f, dt, y)
    ys(:, i+1) = y(:)
end do

end subroutine

function f(x) result(x_new)
! Callback function with derivatives of the Lorentz equations
real(dp), parameter :: s = 10.0_dp, b = 8.0_dp/3, r = 28.0_dp
real(dp), intent(in) :: x(3)
real(dp) :: x_new(3)

x_new(1) = s*(x(2) - x(1))
x_new(2) = r*x(1) - x(2) - x(1)*x(3)
x_new(3) = x(1)*x(2) - b*x(3)
end function

end module