module lorentz
use types, only: dp
use integrators, only: euler
implicit none

contains

subroutine integrate(ys)
! Integration parameters
real(dp), parameter :: t_start = 0.0_dp, t_stop = 1.0_dp, dt = 0.01_dp
integer, parameter :: steps = 1000
integer :: i
real(dp), intent(out) :: ys(3, steps)

! Initial conditions
real(dp) :: y(3)
y = [13.0_dp, 8.1_dp, 45.0_dp]

! Integrator loop
do i = 0, steps-1
    y = euler(f, dt, y, t_start)
    ys(:, i+1) = y(:)
end do

end subroutine

function f(x) result(x_new)
! Callback function for the Lorentz equations
real(dp), parameter :: s = 10.0_dp, b = 8.0_dp/3, r = 28.0_dp
real(dp), intent(in) :: x(3)
real(dp) :: x_new(3)

x_new(1) = s*(x(2) - x(1))
x_new(2) = r*x(1) - x(2) - x(1)*x(3)
x_new(3) = x(1)*x(2) - b*x(3)
end function

end module