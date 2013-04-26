module integrators
use types, only: dp
implicit none
private
public euler, midpoint, rk4, heun, mod_euler

contains

function euler(f, h, y) result(y_new)
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
k2 = h*f(y + 0.5_dp*k1)
y_new = y + k2
end function

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
k2 = h*f(y + 0.5_dp*k1)
k3 = h*f(y + 0.5_dp*k2)
k4 = h*f(y + k3)
y_new = y + k1/6.0_dp + k2/3.0_dp + k3/3.0_dp + k4/6.0_dp
end function

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
k2 = h*f(y + 2.0_dp*k1/3)
y_new = y + 0.25_dp*k1 + 0.75_dp*k2
end function

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
k2 = h*f(y + k1)
y_new = y + 0.5_dp*k1 + 0.5_dp*k2
end function

end module

! function y=rkckODEinter41(coeff,yold,told)
! %
! % initialize the Cash-Karp coefficients
! % defined in the tableau in lab 4,
! % section 3.5
! %
!   a = [.2, 0.3, 0.6, 1.0, 0.875];
!   c1 = [37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0];
!   c2= [2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, .25];
!   b=ones(5);
!   c2 = c1 - c2;
!   b(1,1) =0.2; 
!   b(2,1)= 3.0/40.0; 
!   b(2,2)=9.0/40.0;
!   b(3,1)=0.3; 
!   b(3,2)=-0.9; 
!   b(3,3)=1.2;
!   b(4,1)=-11.0/54.0; 
!   b(4,2)=2.5; 
!   b(4,3)=-70.0/27.0; 
!   b(4,4)=35.0/27.0;
!   b(5,1)=1631.0/55296.0; 
!   b(5,2)=175.0/512.0; 
!   b(5,3)=575.0/13824.0;
!   b(5,4)=44275.0/110592.0; 
!   b(5,5)=253.0/4096.0;
! 
! % set up arrays
!   
!   derivArray=zeros(6,length(yold));
!   ynext=zeros([1,length(yold)]);
!   bsum=zeros([1,length(yold)]);
!   derivArray(1,:)=derivsinter41(coeff,yold,told);
!   
!   % calculate step
!   
!   y=yold;
!   for i=1:5
!     bsum=0.;
!     for j=1:i
!       bsum=bsum + b(i,j)*derivArray(j,:);
!     end
!     derivArray(i+1,:)=derivsinter41(coeff,y + coeff.dt*bsum,told + a(i)*coeff.dt);
!     ynext = ynext + c1(i)*derivArray(i,:);
!   end
!   y = y + coeff.dt*(ynext + c1(6)*derivArray(6,:));
! 
