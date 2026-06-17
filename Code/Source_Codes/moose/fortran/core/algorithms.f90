!===============================================================================
! A few basic algorithms
!===============================================================================
module moose_algorithms
  use iso_fortran_env
  implicit none



  interface swap
     procedure :: swap_i, swap_r
  end interface swap


  interface binary_search_L
     procedure :: binary_search_L_int
     procedure :: binary_search_L_real64
  end interface
  interface binary_search_R
     procedure :: binary_search_R_int
     procedure :: binary_search_R_real64
  end interface


  interface quicksort
     procedure :: quicksort_array1d
     procedure :: quicksort_array2d
  end interface quicksort


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  elemental subroutine swap_i(a, b)
  integer, intent(inout) :: a, b

  integer :: tmp


  tmp = a;   a = b;   b = tmp

  end subroutine swap_i
  !-----------------------------------------------------------------------------
  elemental subroutine swap_r(a, b)
  real(real64), intent(inout) :: a, b

  real(real64) :: tmp


  tmp = a;   a = b;   b = tmp

  end subroutine swap_r
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function binary_search_L_int(A, T) result(m)
  !
  ! return leftmost/smallest index m with A(m) < T <= A(m+1), in particular
  ! T <= A(1) for m = 0 and A(n) < T for m = n = size(A)
  !
  ! note: m is the rank of T in A, i.e. the number of elements that are < T
  !
  integer, intent(in) :: A(:)
  integer, intent(in) :: T
  integer             :: m

  integer :: l, r


  l = 1
  r = ubound(A,1) + 1
  do while (l < r)
     m = (l+r)/2
     if (A(m) < T) then
        l = m + 1
     else
        r = m
     endif
  enddo
  m = l - 1

  end function binary_search_L_int
  !-----------------------------------------------------------------------------
  pure function binary_search_L_real64(A, T) result(m)
  !
  ! return leftmost/smallest index m with A(m) < T <= A(m+1), in particular
  ! T <= A(1) for m = 0 and A(n) < T for m = n = size(A)
  !
  ! note: m is the rank of T in A, i.e. the number of elements that are < T
  !
  real(real64), intent(in) :: A(:)
  real(real64), intent(in) :: T
  integer                  :: m

  integer :: l, r


  l = 1
  r = ubound(A,1) + 1
  do while (l < r)
     m = (l+r)/2
     if (A(m) < T) then
        l = m + 1
     else
        r = m
     endif
  enddo
  m = l - 1

  end function binary_search_L_real64
  !-----------------------------------------------------------------------------
  pure function binary_search_R_int(A, T) result(m)
  !
  ! return rightmost/largest index m with A(m) <= T < A(m+1), in particular
  ! T < A(1) for m = 0 and A(n) <= T for m = n = size(A)
  !
  integer, intent(in) :: A(:)
  integer, intent(in) :: T
  integer             :: m

  integer :: l, r


  l = 1
  r = ubound(A,1) + 1
  do while (l < r)
     m = (l+r)/2
     if (A(m) > T) then
        r = m
     else
        l = m + 1
     endif
  enddo
  m = l - 1

  end function binary_search_R_int
  !-----------------------------------------------------------------------------
  pure function binary_search_R_real64(A, T) result(m)
  !
  ! return rightmost/largest index m with A(m) <= T < A(m+1), in particular
  ! T < A(1) for m = 0 and A(n) <= T for m = n = size(A)
  !
  real(real64), intent(in) :: A(:)
  real(real64), intent(in) :: T
  integer                  :: m

  integer :: l, r


  l = 1
  r = ubound(A,1) + 1
  do while (l < r)
     m = (l+r)/2
     if (A(m) > T) then
        r = m
     else
        l = m + 1
     endif
  enddo
  m = l - 1

  end function binary_search_R_real64
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function argsort(A)
  !
  ! return indices that would sort the array *A*
  !
  real(real64), intent(in) :: A(:)
  integer                  :: argsort(size(A))

  real(real64), allocatable :: tmp(:)
  integer :: i


  allocate (tmp, source = A)
  argsort = [(i, i=1,size(A))]
  call quicksort_array1d_main(tmp, 1, size(A), argsort)

  end function argsort
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine quicksort_array1d(A, i1, i2)
  real(real64), intent(inout) :: A(:)
  integer,      intent(in   ) :: i1, i2

  integer, allocatable :: sort_order(:)
  integer :: i


  sort_order = [(i, i=1,size(A))]
  call quicksort_array1d_main(A, i1, i2, sort_order)

  end subroutine quicksort_array1d
  !-----------------------------------------------------------------------------
  recursive subroutine quicksort_array1d_main(A, i1, i2, sort_order)
  real(real64), intent(inout) :: A(:)
  integer,      intent(in   ) :: i1, i2
  integer,      intent(inout) :: sort_order(size(A))

  integer, allocatable :: order(:)
  integer :: ip


  if (i1 >= i2) return
  ip = partition(A, i1, i2)
  call quicksort_array1d_main(A, i1, ip - 1, sort_order)
  call quicksort_array1d_main(A, ip, i2, sort_order)

  contains
  !.............................................................................
  function partition(A, i1, i2) result(ip)
  implicit none

  real(real64), intent(inout) :: A(:)
  integer,      intent(in)    :: i1, i2
  integer                     :: ip

  real(real64) :: x
  integer :: i, j


  x = A(i1)
  i = i1 - 1
  j = i2 + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     enddo
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     enddo
     if (i < j) then
        call swap(A(i), A(j))
        call swap(sort_order(i), sort_order(j))
     elseif (i == j) then
        ip = i + 1
        return
     else
        ip = i
        return
     endif
  enddo

  end function partition
  !.............................................................................
  end subroutine quicksort_array1d_main
  !=============================================================================
  recursive subroutine quicksort_array2d(A, i1, i2, jpivot, axis)
  !
  ! sort rows (axis = 1) or columns (axis = 2) of matrix A in range [i1, i2]
  ! based on values in jpivot-th column or row, respectively
  !
  use moose_error, only: VALUE_ERROR
  real(real64), intent(inout) :: A(:,:)
  integer,      intent(in   ) :: i1, i2, jpivot, axis

  integer :: ip


  if (axis < 1  .or.  axis > 2) call VALUE_ERROR("invalid axis", "quicksort")
  if (jpivot < 1  .or.  jpivot > size(A,3-axis)) call VALUE_ERROR("jpivot out of range", "quicksort")
  if (i1 < i2) then
     ip = partition(A, i1, i2, jpivot)
     call quicksort(A, i1, ip - 1, jpivot, axis)
     call quicksort(A, ip, i2,     jpivot, axis)
  endif

  contains
  !.............................................................................
  pure function element(i, j)
  integer, intent(in) :: i, j
  real(real64)        :: element


  if (axis == 1) then
     element = A(i,j)
  else
     element = A(j,i)
  endif

  end function element
  !.............................................................................
  subroutine swap(i, j)
  integer, intent(in) :: i, j

  real(real64) :: tmp(size(A,3-axis)), x


  if (axis == 1) then
     tmp    = A(i,:)
     A(i,:) = A(j,:)
     A(j,:) = tmp
  else
     tmp    = A(:,i)
     A(:,i) = A(:,j)
     A(:,j) = tmp
  endif

  end subroutine swap
  !.............................................................................
  function partition(A, i1, i2, jpivot) result(ip)
  implicit none

  real(real64), intent(inout) :: A(:,:)
  integer,      intent(in)    :: i1, i2, jpivot
  integer                     :: ip

  real(real64) :: x
  integer :: i, j


  x = element(i1, jpivot)
  i = i1 - 1
  j = i2 + 1

  do
     j = j-1
     do
        if (element(j, jpivot) <= x) exit
        j = j-1
     enddo
     i = i+1
     do
        if (element(i, jpivot) >= x) exit
        i = i+1
     enddo
     if (i < j) then
        call swap(i, j)
     elseif (i == j) then
        ip = i + 1
        return
     else
        ip = i
        return
     endif
  enddo

  end function partition
  !.............................................................................
  end subroutine quicksort_array2d
  !=============================================================================


  !=============================================================================
  ! Intersect L1 (line through p1 in direction of v1) with
  !           L2 (line through p2 in direction of v2)
  ! Return coordinates s1 and s2 along L1 and L2, respectively, and the
  ! intersection point x in real space.
  ! In case lines are parallel: ierr = 1
  ! if in addition lines are identical: set intersection point to p2
  !=============================================================================
  subroutine intersect_lines_2D(p1, v1, p2, v2, s1, s2, x, ierr)
  real(real64), intent(in)  :: p1(2), v1(2), p2(2), v2(2)
  real(real64), intent(out) :: s1, s2, x(2)
  integer,      intent(out) :: ierr

  real(real64) :: det, dp(2), u(2), p, d


  ierr = 0
  det  = -v1(1)*v2(2) + v1(2)*v2(1)
  dp   = p2 - p1

  ! check if lines are parallel
  p    = max(maxval(p1),maxval(p2))
  if (abs(det) < epsilon(p)) then
     u(1) = -v1(2)
     u(2) =  v1(1)
     d    = sum(u * dp)

     ! lines are the same
     if (abs(d) < epsilon(p)) then
        x  = p2
        s2 = 0.d0
        s1 = sum(dp*v1) / sqrt(sum(v1**2))
     endif

     ierr = 1
     return
  endif

  ! lines are not parallel
  s1 = (-v2(2)*dp(1) + v2(1)*dp(2)) / det
  s2 = (-v1(2)*dp(1) + v1(1)*dp(2)) / det
  x  = p1 + s1 * v1

  end subroutine intersect_lines_2D
  !=============================================================================


  !=============================================================================
  function xsegments(p1, q1, p2, q2)
  !
  ! check if segment p1 - >q1 intersects with segment p2 -> q2
  !
  real(real64), intent(in)  :: p1(2), q1(2), p2(2), q2(2)
  logical                   :: xsegments

  real(real64) :: s1, s2, x(2)
  integer :: ipara


  xsegments = .false.
  call intersect_lines_2D(p1, q1-p1, p2, q2-p2, s1, s2, x, ipara)
  if (ipara == 1) return
  if (s1 < 0.d0  .or.  s1 > 1.d0) return
  if (s2 < 0.d0  .or.  s2 > 1.d0) return
  xsegments = .true.

  end function xsegments
  !=============================================================================


  !=============================================================================
  pure subroutine solve_quadratic_equation(a, b, c, xmin, xmax, x, n)
  !
  ! calculate roots of a*x^2 + b*x + c, with boundary condition xmin <= x <= xmax
  !
  real(real64), intent(in)  :: a, b, c, xmin, xmax
  real(real64), intent(out) :: x(2)
  integer,      intent(out) :: n

  real(real64), parameter :: eps = 1.d-8

  real(real64) :: b4ac, D, bsgna, mod2a, modb, modc, xA, xB, xmin2a, xmax2a, x2


  n = 0
  ! 1. there are no roots ..............................................
  b4ac  = b**2 - 4.d0 * a * c
  if (b4ac < 0.d0) return


  ! 2. c = 0 ...........................................................
  if (c == 0.d0) then
     ! x1 = 0 is a root, but is it in [xmin, xmax]?
     x = 0.d0
     if (xmin <= x(1)  .and.  x(1) <= xmax) n = 1

     ! x2 = -b/a is another root, but only if a /= 0
     if (a /= 0.d0) then
        x2 = -b / a
        ! check if x2 is in [xmin, xmax]
        if (xmin <= x2  .and.  x2 <= xmax) then
           if (0.d0 < x(1)  .and.  n == 1) then
              x(2) = x2
           else
              x(1) = x2
           endif
           n = n + 1
        endif
     endif
     return
  endif
  ! from now on: c /= 0 ................................................


  ! 3. |a|, |b| << |c| (xmin, xmax should be of order 1) ...............
  mod2a = 2.d0 * abs(a)
  modb  = abs(b)
  modc  = abs(c)
  if (mod2a < eps * modc  .and.  modb < eps * modc) return


  ! 4. |a| << |b| (including a = 0, but not a = b = 0) .................
  if (mod2a < eps * modb) then
     if (modc > modb) return

     x(1) = -c/b * (1.d0 + a/b/4)
     if (xmin <= x(1)  .and.  x(1) <= xmax) n = 1
     return
  endif


  ! 5. evaluate roots ..................................................
  D      = sqrt(b4ac)
  bsgna  = b;   if (a >= 0.d0) bsgna = -b
  xmin2a = mod2a * xmin
  xmax2a = mod2a * xmax
  xA     = bsgna + D
  xB     = bsgna - D
  ! no roots in [xmin, xmax]
  if (xA < xmin2a) return
  if (xB > xmax2a) return

  ! smaller root < xmin
  if (xB < xmin2a) then
     ! larger root > xmax
     if (xA > xmax2a) return

     ! larger root in [xmin, xmax]
     x(1) = xA / mod2a
     n    = 1

  ! smaller root in [xmin, xmax]
  else
     x(1) = xB / mod2a
     n    = 1

     ! larger root also in [xmin, xmax]
     if (xA <= xmax2a) then
        x(2) = xA / mod2a
        n    = 2
     endif
  endif

  end subroutine solve_quadratic_equation
  !=============================================================================


  !-----------------------------------------------------------------------------
  pure function fsal(x, dim)
  !
  ! check if first vector element along dimension *dim* is same as last
  !
  real(real64), intent(in) :: x(:,:)
  integer,      intent(in) :: dim
  logical                  :: fsal


  fsal = .false.
  if (dim == 1) then
     fsal = all(abs(x(:,lbound(x,2)) - x(:,ubound(x,2))) < 1.d2 * maxval(abs(x)) * epsilon(1.d0))
  elseif (dim == 2) then
     fsal = all(abs(x(lbound(x,1),:) - x(ubound(x,1),:)) < 1.d2 * maxval(abs(x)) * epsilon(1.d0))
  endif

  end function fsal
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function reverse_array(A) result(Aout)
  real(real64), intent(in) :: A(:)
  real(real64)             :: Aout(lbound(A,1):ubound(A,1))

  integer :: i, i1, i2


  i1 = lbound(A,1)
  i2 = ubound(A,1)
  do i=i1,i2
     Aout(i2+i1-i) = A(i)
  enddo

  end function reverse_array
  !-----------------------------------------------------------------------------
  pure function reverse_matrix(M, axis) result(Mout)
  real(real64), intent(in) :: M(:,:)
  integer,      intent(in) :: axis
  real(real64)             :: Mout(lbound(M,1):ubound(M,1), lbound(M,2):ubound(M,2))


  if (axis == 1) then
     Mout = reverse_rows(M)
  elseif (axis == 2) then
     Mout = reverse_columns(M)
  endif

  end function reverse_matrix
  !-----------------------------------------------------------------------------
  pure function reverse_rows(M) result(Mout)
  real(real64), intent(in) :: M(:,:)
  real(real64)             :: Mout(lbound(M,1):ubound(M,1), lbound(M,2):ubound(M,2))

  integer :: i, i1, i2


  i1 = lbound(M,1)
  i2 = ubound(M,1)
  do i=i1,i2
     Mout(i2+i1-i,:) = M(i,:)
  enddo

  end function reverse_rows
  !-----------------------------------------------------------------------------
  pure function reverse_columns(M) result(Mout)
  real(real64), intent(in) :: M(:,:)
  real(real64)             :: Mout(lbound(M,1):ubound(M,1), lbound(M,2):ubound(M,2))

  integer :: j, j1, j2


  j1 = lbound(M,2)
  j2 = ubound(M,2)
  do j=j1,j2
     Mout(:,j2+j1-j) = M(:,j)
  enddo

  end function reverse_columns
  !-----------------------------------------------------------------------------

end module moose_algorithms
