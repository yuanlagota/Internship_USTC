module moose_linalg
  use iso_fortran_env
  implicit none


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine inverse_2d(M, M1, istat)
  !
  ! calculate inverse of M, or return istat = 1 if it does not exist
  !
  real(real64), intent(in)  :: M(:,:)
  real(real64), intent(out) :: M1(size(M,1), size(M,2))
  integer,      intent(out) :: istat

  real(real64) :: Mdet


  istat = 1
  if (size(M,1) == 2  .and.  size(M,2) == 2) then
     Mdet = M(1,1) * M(2,2)  -  M(1,2) * M(2,1)
     if (Mdet == 0.d0) return

     M1(1,1) =  M(2,2) / Mdet
     M1(2,2) =  M(1,1) / Mdet
     M1(1,2) = -M(1,2) / Mdet
     M1(2,1) = -M(2,1) / Mdet
     istat   = 0
  endif

  end subroutine inverse_2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function second_partial_derivative_test_2d(H) result(icase)
  !
  ! perform test for Hessian matrix H, D = det(H): return icase
  !    = 1, if this is a local minimum
  !    = 2, if this is a local maximum 
  !    = 3, if this is saddle point
  !    = 4, if test is inconclusive
  !
  real(real64), intent(in) :: H(2,2)
  integer                  :: icase

  real(real64) :: D


  D = H(1,1)*H(2,2) - H(1,2)*H(2,1)
  if (D > 0.d0  .and.  H(1,1) > 0.d0) then
     icase = 1
  elseif (D > 0.d0  .and.  H(1,1) < 0.d0) then
     icase = 2
  elseif (D < 0.d0) then
     icase = 3
  else
     icase = 4
  endif

  end function second_partial_derivative_test_2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine hessian2d_analysis(Hxx, Hxy, Hyy, l1, l2, v1, v2, alpha)
  !
  ! compute eigenvalues l1,l2 and eigenvectors v1,v2 of 2-D Hessian
  !
  ! v1 will point in positive x-direction (or positive [cos(alpha), sin(alpha)] direction),
  ! and v2 will have a counter-clockwise orientation
  !
  real(real64), intent(in   ) :: Hxx, Hxy, Hyy
  real(real64), intent(  out) :: l1, l2, v1(2), v2(2)
  real(real64), intent(in   ), optional :: alpha

  real(real64) :: ac2, ac4, b2, e1(2), e2(2)


  ! compute eigenvalues
  ac2 = 0.5d0  * (Hxx + Hyy)
  ac4 = 0.25d0 * (Hxx - Hyy)**2
  b2  = Hxy**2
  l1  = ac2 - dsqrt(ac4 + b2)
  l2  = ac2 + dsqrt(ac4 + b2)


  e1 = [1.d0, 0.d0];   if (present(alpha)) e1 = [cos(alpha), sin(alpha)]
  e2 = [-e1(2), e1(1)]
  ! construct eigenvectors
  v1 = make_eigenvector(l1, e1)
  v2 = make_eigenvector(l2, e2)

  contains
  !.............................................................................
  function make_eigenvector(l, e) result(v)
  real(real64), intent(in) :: l, e(2)
  real(real64)             :: v(2)

  real(real64) :: a1, a2, a3, aa1, aa2, aa3


  a1 = Hxx - l;   aa1 = abs(a1)
  a2 = Hxy;       aa2 = abs(a2)
  a3 = Hyy - l;   aa3 = abs(a3)
  if (aa1 > max(aa2,aa3)) then
     v = [-a2/a1, 1.d0]
  elseif (aa3 > max(aa1,aa2)) then
     v = [1.d0, -a2/a3]
  elseif (aa1 > aa3) then
     v = [-a3/a2, 1.d0]
  else
     v = [1.d0, -a1/a2]
  endif

  ! normalize
  v = v / sqrt(sum(v**2))

  ! adjust oriantation
  if (sum(v*e) < 0.d0) v = - v

  end function make_eigenvector
  !.............................................................................
  end subroutine hessian2d_analysis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine stability_analysis(H, lambda1, lambda2, v1, v2, theta)
  !
  ! calculate eigenvalues and eigenvectors for stable and unstable manifolds
  ! NOTE: eigenvectors are oriented in positive x-direction (or theta-direction, if present)
  !
  ! lambda1, v1: unstable direction
  ! lambda2, v2: stable direction
  !
  use moose_math, only: pi
  real(real64), intent(in   ) :: H(2,2)
  real(real64), intent(  out) :: lambda1, lambda2, v1(2), v2(2)
  real(real64), intent(in   ), optional :: theta


  real(real64) :: A(2,2), Q, phi(2), l(2), delta
  integer      :: i


  A(1,1) = -H(1,2);   A(2,1) = H(1,1)
  A(1,2) = -H(2,2);   A(2,2) = H(2,1)

  Q = - (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  if (Q < 0.d0) then
     return
  endif
  ! calculate eigenvalues
  lambda1 =  sqrt(Q);  l(1) = lambda1
  lambda2 = -sqrt(Q);  l(2) = lambda2

  ! calculate eigenvectors
  do i=1,2
     phi(i) = atan2(-A(1,1) + l(i), A(1,2))
     delta  = phi(i);   if (present(theta)) delta = phi(i) - theta
     if (abs(delta) > pi/2.d0  .and.  abs(delta) < 3.d0*pi/2.d0) then
        phi(i) = phi(i) - sign(pi, delta)
     endif
  enddo
  v1(1) = cos(phi(1))
  v1(2) = sin(phi(1))
  v2(1) = cos(phi(2))
  v2(2) = sin(phi(2))

  end subroutine stability_analysis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine linregress(x, y, a, b)
  real(real64), intent(in   ) :: x(:), y(size(x))
  real(real64), intent(  out) :: a, b

  real(real64) :: xsum, x2sum, xysum, ysum
  integer :: n


  n = size(x)
  xsum = sum(x)
  x2sum = sum(x**2)
  xysum = sum(x*y)
  ysum = sum(y)

  a = (n * xysum - xsum * ysum) / (n * x2sum - xsum**2)
  b = (ysum - a * xsum) / n

  end subroutine linregress
  !-----------------------------------------------------------------------------

end module moose_linalg
