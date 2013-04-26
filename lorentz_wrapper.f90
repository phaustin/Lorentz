module lorentz_wrapper
use lorentz, only : butterfly

implicit none

contains

subroutine c_butterfly() bind (c)
call butterfly()
end subroutine

end module
