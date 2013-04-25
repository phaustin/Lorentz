module types
implicit none
private
public sp, dp

integer, parameter :: sp = kind(0.0)      ! single precision
integer, parameter :: dp = kind(0.d0)     ! double precision

end module