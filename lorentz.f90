program lorentz
use types, only: dp
use utils
implicit none

! read a file and print its contents
integer :: u
real(dp) :: a(1116)
real(dp) :: B(279, 4)

open(newunit=u, file="data.txt", form="formatted", status="old", action="read")
read(u, *) a
B = reshape(a, [279, 4], order=[2, 1])
close(u)
    
end program