! Runge-Kutta method Tableaus. Recall Fortran is column major.
module tableaus
use types, only: dp
implicit none
private
public a_mp, c_mp, b_mp
public a_hm, c_hm, b_hm
public a_me, c_me, b_me
public a_rk4, c_rk4, b_rk4
public a_ck, c1_ck, c2_ck, b_ck

! Midpoint method tableau
real(dp), parameter :: a_mp(2) = [0.0_dp, 0.5_dp], c_mp(2) = [0.0_dp, 1.0_dp]
real(dp), parameter :: b_mp(2, 2) = &
    reshape([0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp], [2, 2])

! Heun's method tableau
real(dp), parameter :: a_hm(2) = [0.0_dp, 2.0_dp/3] 
real(dp), parameter :: c_hm(2) = [1.0_dp/4, 3.0_dp/4]
real(dp), parameter :: b_hm(2, 2) = &
    reshape([0.0_dp, 2.0_dp/3, 0.0_dp, 0.0_dp], [2, 2])

! Modified Euler method tableau
real(dp), parameter :: a_me(2) = [0.0_dp, 1.0_dp] 
real(dp), parameter :: c_me(2) = [0.5_dp, 0.5_dp]
real(dp), parameter :: b_me(2, 2) = &
    reshape([0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], [2, 2])

! Fourth-order Runge-Kutta method tableau
real(dp), parameter :: a_rk4(4) = [0.0_dp, 0.5_dp, 0.5_dp, 1.0_dp]
real(dp), parameter :: c_rk4(4) = [1.0_dp/6, 1.0_dp/3, 1.0_dp/3, 1.0_dp/6]
real(dp), parameter :: b_rk4(4, 4) = reshape([ &
    0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
    [4, 4])

! Cash-Karp method (c1 is fifth order and c2 is fourth order)
real(dp), parameter :: a_ck(6) = [0.0_dp, 0.2_dp, 0.3_dp, 0.6_dp, 1.0_dp, &
    0.875_dp]
real(dp), parameter :: c1_ck(6) = [37.0_dp/378, 0.0_dp, 250.0_dp/621, &
    125.0_dp/594, 0.0_dp, 512.0_dp/1771]
real(dp), parameter :: c2_ck(6) = [2825.0_dp/27648, 0.0_dp, &
    18575.0_dp/48384, 13525.0_dp/55296, 277.0_dp/14336, 0.25_dp]
real(dp), parameter :: b_ck(6, 6) = reshape([ &
    0.0_dp, 0.2_dp, 3.0_dp/40, 3.0_dp/10, -11.0_dp/54, 1631.0_dp/55296, &
    0.0_dp, 0.0_dp, 9.0_dp/40, -0.9_dp, 2.5_dp, 175.0_dp/512, &
    0.0_dp, 0.0_dp, 0.0_dp, 6.0_dp/5, -70.0_dp/27, 575.0_dp/13824, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 35.0_dp/27, 44275.0_dp/110592, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 253.0_dp/4096, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
    [6, 6])

end module