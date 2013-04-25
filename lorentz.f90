module lorentz
use types, only: dp
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

end module