module lorentz
use types, only: dp
use integrators, only: euler
implicit none

contains

subroutine assimilate_data(B)
! read a file and make its contents available to Python
integer :: u
real :: a(1116)
real, intent(out) :: B(279, 4)

open(newunit=u, file="data.txt", form="formatted", status="old", action="read")
read(u, *) a
B = reshape(a, [279, 4], order=[2, 1])
close(u)
end subroutine

subroutine butterfly()
! Integration parameters
real(dp), parameter :: t_start = 0.0_dp, t_stop = 0.1_dp, dt = 0.01_dp

! Initial conditions
real(dp) :: y(3)
y = [13.0_dp, 8.1_dp, 45.0_dp]

!Integrator loop; loops are overrated, just STEP; frikking woot.
y = euler(f, dt, y, t_start)
print*, y
y = euler(f, dt, y, t_start+dt)
print*, y
y = euler(f, dt, y, t_start+2*dt)
print*, y
y = euler(f, dt, y, t_start+3*dt)
print*, y

end subroutine

function f(x) result(x_new)
! Parameters for the Lorentz equations
real(dp), parameter :: s = 10.0_dp, b = 8.0_dp/3, r = 28.0_dp
real(dp), intent(in) :: x(3)
real(dp) :: x_new(3)

x_new(1) = s*(x(2) - x(1))
x_new(2) = r*x(1) - x(2) - x(1)*x(3)
x_new(3) = x(1)*x(2) - b*x(3)
end function

end module