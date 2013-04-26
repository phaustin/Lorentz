! Library of integrators that use callbacks to integrate aribitrary functions.
! All methods are explicit and their tableaus are assumed to be lower 
! triangular to simplify coding.
module integrators
use types, only: dp
use tableaus, only: a_mp, c_mp, b_mp
use tableaus, only: a_rk4, c_rk4, b_rk4
use tableaus, only: a_hm, c_hm, b_hm
use tableaus, only: a_me, c_me, b_me
use tableaus, only: a_ck, c1_ck, c2_ck, b_ck
implicit none
private
public fwd_euler, midpoint, rk4, heun, mod_euler, ck4, ck5

contains

! Forward Euler method
function fwd_euler(f, h, y) result(y_new)
real(dp), intent(in) :: h
real(dp), intent(in) :: y(3)
real(dp) :: y_new(3)
interface
    function f(x) result(x_new)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(3)
    real(dp) :: x_new(3)
    end function
end interface
y_new = y + h*f(y)
end function

! Midpoint method
function midpoint(f, h, y) result(y_new)
real(dp), intent(in) :: h
real(dp), intent(in) :: y(3)
real(dp) :: y_new(3), k1(3), k2(3)
interface
    function f(x) result(x_new)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(3)
    real(dp) :: x_new(3)
    end function
end interface
k1 = h*f(y)
k2 = h*f(y + b_mp(2, 1)*k1)
y_new = y + c_mp(1)*k1 + c_mp(2)*k2
end function

! Heun's method
function heun(f, h, y) result(y_new)
real(dp), intent(in) :: h
real(dp), intent(in) :: y(3)
real(dp) :: y_new(3), k1(3), k2(3)
interface
    function f(x) result(x_new)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(3)
    real(dp) :: x_new(3)
    end function
end interface
k1 = h*f(y)
k2 = h*f(y + b_mp(2, 1)*k1)
y_new = y + c_hm(1)*k1 + c_hm(2)*k2
end function

! Modified Euler's method
function mod_euler(f, h, y) result(y_new)
real(dp), intent(in) :: h
real(dp), intent(in) :: y(3)
real(dp) :: y_new(3), k1(3), k2(3)
interface
    function f(x) result(x_new)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(3)
    real(dp) :: x_new(3)
    end function
end interface
k1 = h*f(y)
k2 = h*f(y + b_me(2, 1)*k1)
y_new = y + c_me(1)*k1 + c_me(2)*k2
end function

! Fourth-order Runge-Kutta method
function rk4(f, h, y) result(y_new)
real(dp), intent(in) :: h
real(dp), intent(in) :: y(3)
real(dp) :: y_new(3), k1(3), k2(3), k3(3), k4(3)
interface
    function f(x) result(x_new)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(3)
    real(dp) :: x_new(3)
    end function
end interface
k1 = h*f(y)
k2 = h*f(y + b_rk4(2, 1)*k1)
k3 = h*f(y + b_rk4(3, 1)*k1 + b_rk4(3, 2)*k2)
k4 = h*f(y + b_rk4(4, 1)*k1 + b_rk4(4, 2)*k2 + b_rk4(4, 3)*k3)
y_new = y + c_rk4(1)*k1 + c_rk4(2)*k2 + c_rk4(3)*k3 + c_rk4(4)*k4
end function

! Fourth-order Cash-Karp method
function ck4(f, h, y) result(y_new)
real(dp), intent(in) :: h
real(dp), intent(in) :: y(3)
real(dp) :: y_new(3), k1(3), k2(3), k3(3), k4(3), k5(3), k6(3)
interface
    function f(x) result(x_new)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(3)
    real(dp) :: x_new(3)
    end function
end interface
k1 = h*f(y)
k2 = h*f(y + b_ck(2, 1)*k1)
k3 = h*f(y + b_ck(3, 1)*k1 + b_ck(3, 2)*k2)
k4 = h*f(y + b_ck(4, 1)*k1 + b_ck(4, 2)*k2 + b_ck(4, 3)*k3)
k5 = h*f(y + b_ck(5, 1)*k1 + b_ck(5, 2)*k2 + b_ck(5, 3)*k3 + b_ck(5, 4)*k4)
k6 = h*f(y + b_ck(6, 1)*k1 + b_ck(6, 2)*k2 + b_ck(6, 3)*k3 + b_ck(6, 4)*k4 + &
    b_ck(6, 5)*k5)
y_new = y + c2_ck(1)*k1 + c2_ck(2)*k2 + c2_ck(3)*k3 + c2_ck(4)*k4 + &
    c2_ck(6)*k4 + c2_ck(6)*k6
end function

! Fifth-order Cash-Karp method
function ck5(f, h, y) result(y_new)
real(dp), intent(in) :: h
real(dp), intent(in) :: y(3)
real(dp) :: y_new(3), k1(3), k2(3), k3(3), k4(3), k5(3), k6(3)
interface
    function f(x) result(x_new)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: x(3)
    real(dp) :: x_new(3)
    end function
end interface
k1 = h*f(y)
k2 = h*f(y + b_ck(2, 1)*k1)
k3 = h*f(y + b_ck(3, 1)*k1 + b_ck(3, 2)*k2)
k4 = h*f(y + b_ck(4, 1)*k1 + b_ck(4, 2)*k2 + b_ck(4, 3)*k3)
k5 = h*f(y + b_ck(5, 1)*k1 + b_ck(5, 2)*k2 + b_ck(5, 3)*k3 + b_ck(5, 4)*k4)
k6 = h*f(y + b_ck(6, 1)*k1 + b_ck(6, 2)*k2 + b_ck(6, 3)*k3 + b_ck(6, 4)*k4 + &
    b_ck(6, 5)*k5)
y_new = y + c1_ck(1)*k1 + c1_ck(2)*k2 + c1_ck(3)*k3 + c1_ck(4)*k4 + &
    c1_ck(6)*k4 + c1_ck(6)*k6
end function

end module
