module moose_contours
  use iso_fortran_env
  use moose_txtio, only: txtio, init_txtio
  use moose_rlist
  use moose_mfunc
  implicit none
  private


  integer, parameter, public :: &
     XCONTOUR_STABLE_DIRECTION     =  1, & ! branch of the stable (incoming) manifold
     XCONTOUR_UNSTABLE_DIRECTION   = -1, & ! branch of the unstable (outgoing) manifold
     XCONTOUR_POSITIVE_ORIENTATION =  1, & ! branch in positive reference direction
     XCONTOUR_NEGATIVE_ORIENTATION = -1, & ! branch in negative reference direction
     XCONTOUR_BRANCH_DIRECTION(0:3)   = [-1, -1, 1, 1], &
     XCONTOUR_BRANCH_ORIENTATION(0:3) = [-1, 1, -1, 1]



  ! contour lines ..............................................................
  type, extends(txtio), public :: contour
     type(rlist), pointer :: L
     real(real64), pointer :: x(:,:), t(:,:), dpsi(:)

     contains
     procedure :: free
     procedure :: npoints, point, add_point
     procedure :: arclengths
     procedure :: interp
     procedure :: write_formatted
  end type contour


  interface contour
     procedure :: generate_contour
  end interface
  ! contour ....................................................................



  ! contour lines from hyperbolic point ........................................
  type, public :: xcontour
     type(contour) :: branch(0:3)
     real(real64) :: x(2)
     integer      :: iconnect(0:3)
     real(real64) :: uconnect(0:3)

     contains
     procedure :: savetxt
  end type xcontour


  interface xcontour
     procedure :: generate
  end interface
  ! xcontour ...................................................................



  public :: &
     XCONTOUR_BRANCH, interp_contour

  contains
  !-----------------------------------------------------------------------------


! type contour =================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function generate_contour(Phi, x0, dir, dxmax, istat, &
                        Xpoints, boundary, reverse, epsabs, &
                        dxpmax, nmax, alpha, iconnect, uconnect) result(this)
  !
  ! return contour line Phi = const through x0 in direction 'dir' with steps <= dxmax
  !
  ! return status: istat =
  !    0	successful generation of closed contour line
  !    1	required accuracy not achieved in last step after max. number of correction steps 
  !   -1	contour line out of bounds
  !   -2        contour line connects to hyperbolic point
  !   -3        contour line connects to boundary
  !
  ! optional input:
  !    Xpoints	list of hyperbolic points that contour line may connect to
  !    boundary
  !
  ! numerical parameters:
  !    epsabs	required absolute accuracy
  !    dxpmax	max. allowed distance to hyperbolic point before cut-off
  !    nmax	max. number of correction steps
  !    alpha	damping factor applied to correction steps
  !
  ! optional output:
  !    iconnect     id of hyperbolic point that contour line connects to (istat = -2) or
  !                 id of boundary element that contour line connects to (istat = -3)
  !    uconnect     coordinates of boundary intersection point (istat = -3)
  !
  use moose_math, only: r90
  use moose_analysis, only: stability_analysis
  use moose_hypersurface
  class(scalar_mfunc), intent(in)  :: Phi
  real(real64),        intent(in)  :: x0(2), dxmax
  integer,             intent(in)  :: dir
  integer,             intent(out) :: istat
  type(rlist2),        intent(in),  optional :: Xpoints
  type(hypersurf2d),   intent(in),  optional :: boundary
  logical,             intent(in),  optional :: reverse
  real(real64),        intent(in),  optional :: epsabs, dxpmax, alpha
  integer,             intent(in),  optional :: nmax
  integer,             intent(out), optional :: iconnect
  real(real64),        intent(out), optional :: uconnect
  type(contour)                    :: this

  real(real64), parameter :: dtheta = 0.01d0

  real(real64), allocatable :: xp(:,:), v(:,:,:), lambda(:,:)
  real(real64) :: Phi0, Phi1, x1(2), x2(2), dx, N0(2), N(2), T0(2), T(2), H(2,2), THT, Nmod, Nmod0
  real(real64) :: dxp, d1, dx1, xx(2), tt, uu(1), b, a1, eps1, dxv(2)
  logical :: check_boundary, first_step, prepend
  integer :: i, idir, j, nx, nn, m1


  call assert_mfunc2d(Phi, "contour")


  ! set numerical parameters
  eps1 = 1.d-8;   if (present(epsabs))      eps1 = epsabs      ! required absolute accuracy
  m1   = 16;      if (present(nmax))        m1   = nmax        ! max. number of correction steps
  a1   = 1.d0;    if (present(alpha))       a1   = alpha       ! damping factor applied to correction steps
  idir = 1;       if (dir < 0)              idir = -1          ! trace direction
  nx   = 0;       if (present(Xpoints))     nx   = Xpoints%nelements()
  d1   = dxmax/2; if (present(dxpmax))      d1   = dxpmax
  check_boundary = .false.;   if (present(boundary)) check_boundary = .true.
  prepend        = .false.;   if (present(reverse)) prepend = reverse
  ! internal representation of X-points (hyperbolic points)
  if (nx > 0) then
     allocate (xp(2,0:nx-1), lambda(2,0:nx-1), v(2,2,0:nx-1))
     do i=0,nx-1
        xp(:,i) = Xpoints%element(i)
        H = Phi%hessian(xp(:,i))
        call stability_analysis(H, lambda(1,i), lambda(2,i), v(:,1,i), v(:,2,i))
     enddo
  endif


  ! initialize contour line
  call init_txtio(this, "contour")
  Phi0 = Phi%eval(x0)
  N0 = Phi%deriv(x0);   Nmod0 = norm2(N0)
  T0 = r90(N0) / Nmod0
  allocate (this%L, source=rlist(5))
  call this%add_point(x0, idir*T0, Nmod0, prepend)
  if (present(iconnect)) iconnect = 0
  if (present(uconnect)) uconnect = 0.d0


  ! trace contour line
  x2 = x0
  x1 = x2
  N = N0;   Nmod = Nmod0
  T = T0
  istat = 0
  first_step = .true.
  trace: do
     ! evaluate first and second order derivatives
     H    = Phi%hessian(x2)
     THT  = N(2)**2 * H(1,1) - 2*N(1)*N(2)*H(1,2) + N(1)**2 * H(2,2)


     ! set step size
     b = dtheta * Nmod**3
     if (b < abs(THT) * dxmax) then
        dx = b / abs(THT)
     else
        dx = dxmax
     endif
     ! check distance to X-points
     do j=0,nx-1
        dxp = sqrt(sum((x2-xp(:,j))**2))

        ! contour line connects to X-point?
        if (dxp < d1) then
           if (present(iconnect)) iconnect = j
           xx = xp(:,j) - x1
           dxv(1) = dot_product(xx, v(:,1,j))
           dxv(2) = dot_product(xx, v(:,2,j))
           i = 1;   if (abs(dxv(1)) < abs(dxv(2))) i = 2
           T = sign(1.d0, dxv(i)) * v(:,i,j)
           call this%add_point(xp(:,j), T, 0.d0, prepend)
           istat = -2;   exit trace
        endif

        ! reduce step size near X-point
        if (dxp/2 < dx) dx = dxp/2
     enddo
     if (first_step) then
        ! save first step size
        dx1 = dx
        first_step = .false.
     else
     ! contour connects back to itself?
        dxp = sqrt(sum((x2-x0)**2))
        if (dxp < dx1/4) then
           call this%add_point(x0, idir*T0, Nmod0, prepend)
           exit trace
        endif

        if (dxp/2 < dx) dx = dxp/2
     endif


     ! predictor step
     x2 = x2 + idir * dx * (T - dtheta * N / Nmod / 2)
     ! step out of bounds?
     if (Phi%out_of_bounds(x2)) then
        istat = -1;   exit trace
     endif

     ! corrector step(s)
     do j=1,m1
        Phi1 = Phi%eval(x2)
        N    = Phi%deriv(x2)
        Nmod = norm2(N)
        dx   = (Phi0 - Phi1) / Nmod
        x2   = x2 + a1 * dx * N / Nmod
        ! reached required accuracy?
        if (abs(dx) < eps1) exit
        ! step out of bounds?
        if (Phi%out_of_bounds(x2)) then
           istat = -1;   exit trace
        endif

        ! reached max. number of correction steps
        if (j == m1) then
           istat = 1;   exit trace
        endif
     enddo


     ! check intersection with boundary
     if (check_boundary) then
        if (boundary%intersect(x1, x2, xx, tt, nn, uu)) then
           if (present(iconnect)) iconnect = nn
           if (present(uconnect)) uconnect = uu(1)
           N = Phi%deriv(xx);   Nmod = norm2(N)
           T = r90(N) / Nmod
           call this%add_point(xx, idir*T, Nmod, prepend)
           istat = -3;   exit trace
        endif
     endif


     ! save this point and prepare next step
     N = Phi%deriv(x2);   Nmod = norm2(N)
     T = r90(N) / Nmod
     call this%add_point(x2, idir*T, Nmod, prepend)
     x1 = x2
  enddo trace

  end function generate_contour
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(contour), intent(inout) :: this


  call this%L%free()
  deallocate (this%L)
  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function npoints(this)
  class(contour), intent(in) :: this
  integer                    :: npoints


  npoints = this%L%nelements()

  end function npoints
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function point(this, i) result(x)
  !
  ! pointer to i-th point in contour
  !
  class(contour), intent(in) :: this
  integer,        intent(in) :: i
  real(real64),   pointer    :: x(:)

  real(real64), pointer :: tmp(:)


  tmp => this%L%element(i)
  x => tmp(1:2)

  end function point
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine add_point(this, x, t, dpsi, prepend)
  class(contour), intent(inout) :: this
  real(real64),   intent(in   ) :: x(2), t(2), dpsi
  logical,        intent(in   ) :: prepend


  if (prepend) then
     call this%L%prepend([x, -t, dpsi])
  else
     call this%L%append([x, t, dpsi])
  endif
  this%x => this%L%columns(1, 2)
  this%t => this%L%columns(3, 4)
  this%dpsi => this%L%column(5)

  end subroutine add_point
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function arclengths(this) result(dl)
  class(contour), intent(in) :: this
  real(real64)               :: dl(this%L%nelements()-1)

  integer :: i


  do i=1,this%L%nelements()-1
     ! TODO: approximate arclength based on tangents
     dl(i) = sqrt(sum((this%x(:,i) - this%x(:,i-1))**2))
  enddo

  end function arclengths
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp(this) result(C)
  !
  ! construct interpolated curve from points and tangents along contour
  !
  use moose_math,         only: rowdiv, zero_cumsum
  use moose_interp_curve, only: cubic_hermite_curve, interp_curve
  class(contour), intent(in) :: this
  type(interp_curve)         :: C


  C = cubic_hermite_curve(zero_cumsum(this%arclengths()), transpose(this%x), transpose(this%t))

  end function interp
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(contour),   intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%L%write_formatted(unit, iotype, vlist, iostat, iomsg)

  end subroutine write_formatted
  !---------------------------------------------------------------------
! type contour =================================================================



! type xcontour ================================================================
! auxiliary functions:
  !-----------------------------------------------------------------------------
  function XCONTOUR_BRANCH(direction, orientation) result(k)
  integer, intent(in) :: direction, orientation
  integer             :: k

  integer :: idir


  if (direction == XCONTOUR_STABLE_DIRECTION) then
     idir = 1
  elseif (direction == XCONTOUR_UNSTABLE_DIRECTION) then
     idir = 0
  else
     print 9000, direction
     stop
  endif
 9000 format("ERROR: XCONTOUR_BRANCH called with invalided direction '",i0,"'")


  if (orientation == XCONTOUR_POSITIVE_ORIENTATION) then
     k = 2*idir + 1
  elseif (orientation == XCONTOUR_NEGATIVE_ORIENTATION) then
     k = 2*idir
  else
     print 9001, orientation
     stop
  endif
 9001 format("ERROR: XCONTOUR_BRANCH called with invalid orientation '",i0,"'")

  end function XCONTOUR_BRANCH
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function generate(Psi, Xp, ix, dxmax, theta, boundary, reverse, offset, fX, epsabs, nmax, alpha) result(S)
  !
  ! Construct stable and unstable manifolds for hyperbolic point.
  !
  ! **Parameters:**
  !
  ! :Psi:       Two-dimensional scalar field.
  !
  ! :Xp:        List of (relevant) hyperbolic points in domain of Psi.
  !
  ! :ix:        Index of hyberpolic point in `Xp` for which manifolds should be constructed.
  !
  ! :dxmax:     Max. step size along contour.
  !
  ! **Optional parameters:**
  !
  ! :theta:     Reference direction [rad] for orientation of branches (default: 0).
  !
  ! :boundary:  Contour of domain boundary.
  !
  ! :reverse:   Reverse list of points on invariant manifolds.
  !
  ! :offset:    Initial offset from hyperbolic point for contour tracing (default: `dxmax`/10).
  !
  ! :fX:        Connect to any hyperbolic point within this fraction of `offset` (default: 0.9).
  !
  ! :epsabs:    Required absolute accuracy for points on contour.
  !
  ! :nmax:      Max. number of correction steps for each contour node (default: 16).
  !
  ! :alpha:     Damping factor applied to correction steps (default: 1.0).
  !
  use moose_hypersurface
  use moose_linalg, only: stability_analysis
  class(scalar_mfunc), intent(in)  :: Psi
  type(rlist2),        intent(in)  :: Xp
  integer,             intent(in)  :: ix
  real(real64),        intent(in)  :: dxmax
  real(real64),        intent(in), optional :: theta, offset, fX, epsabs, alpha
  type(hypersurf2d),   intent(in), optional :: boundary
  logical,             intent(in), optional :: reverse
  integer,             intent(in), optional :: nmax
  type(xcontour)                   :: S

  character(len=128) :: filename
  logical      :: reverse_branch
  real(real64) :: lambda1, lambda2, v(2,-1:1), t(2), n(2), xi(2), d0, dxpmax, Psi0, H(2,2), uconnect
  integer      :: istat, iconnect, k, direction, orientation


  call assert_mfunc2d(Psi, "xcontour")
  if (ix < 0  .or.  ix >= Xp%nelements()) then
     write (6, *) "ERROR in xcontour constructor: ix out of range!"
     stop
  endif


  S%x  = Xp%element(ix)
  Psi0 = Psi%eval(S%x)
  H    = Psi%hessian(S%x)
  d0   = 0.1d0 * dxmax;   if (present(offset)) d0     = offset
  dxpmax = 0.9d0 * d0;    if (present(fX))     dxpmax = fX * d0
  call stability_analysis(H, lambda1, lambda2, v(:,-1), v(:,1), theta)


  do direction=-1,1,2
  do orientation=-1,1,2
     k = XCONTOUR_BRANCH(direction, orientation)
     t = orientation * v(:,direction)
     n = [t(2), -t(1)]


     ! 1. move initial point a small step away from X-point
     ! approximate initial point
     xi = S%x + d0 * t
     ! correct initial point to match Psi0
     call root_finder(Psi, xi, Psi0, n, istat, epsabs)
     if (istat /= 0) then
        write (6, *) "ERROR in xcontour: cannot generate initial point!"
        write (6, *) "istat = ", istat
        write (6, *) "PsiT  = ", Psi0
        write (6, *) "Psi1  = ", Psi%eval(S%x + d0 * t)
        write (6, *) "x1    = ", S%x + d0 * t
        stop
     endif


     ! 2. trace separatrix branch
     reverse_branch = direction == 1
     if (present(reverse)) reverse_branch = .not. reverse.eqv.reverse_branch

     S%branch(k) = contour(Psi, xi, -direction, dxmax, istat, &
        Xp, boundary, reverse_branch, epsabs, dxpmax, nmax, alpha, iconnect, uconnect)
     call S%branch(k)%add_point(S%x, -t, 0.d0, .not.reverse_branch)

     select case(istat)
     ! contour line leaves domain of definition
     case(-1)
        S%iconnect(k) = 0

     ! contour line connects to X-point
     case(-2)
        S%iconnect(k) = iconnect + 1

     ! contour line connects to boundary element
     case(-3)
        S%iconnect(k) = -iconnect

     case default
        print 9002, k, istat
        write (filename, 9003) k
        call S%branch(k)%savetxt(filename)
        stop
     end select
     S%uconnect(k) = uconnect
 9002 format("ERROR in xcontour: tracing of branch ",i0," failed with istat = ",i0,"!")
 9003 format("separatrix_branch_",i0,"_ERROR.plt")
  enddo
  enddo

  end function generate
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine savetxt(this, basename, split)
  class(xcontour),  intent(in) :: this
  character(len=*), intent(in) :: basename
  logical,          intent(in), optional :: split

  character(len=len_trim(basename)+7) :: filename
  logical :: split_output
  integer :: k, iu


  split_output = .false.;   if (present(split)) split_output = split


  if (split_output) then
     do k=0,3
        write (filename, 1001) trim(basename), k
        call this%branch(k)%savetxt(filename)
     enddo

  else
     write (filename, 1002) trim(basename)
     open  (newunit=iu, file=filename)
     do k=0,3
        write (iu, '(dt)') this%branch(k)
        write (iu, *)
     enddo
     close (iu)
  endif
 1001 format(a,"_",i0,".txt")
 1002 format(a,".txt")

  end subroutine savetxt
  !-----------------------------------------------------------------------------
! type xcontour ================================================================


! module procedures:
  !-----------------------------------------------------------------------------
  function interp_contour(Phi, x0, dir, dxmax, istat, &
                        Xpoints, boundary, reverse, epsabs, &
                        dxpmax, nmax, alpha, iconnect, uconnect) result(this)
  !
  !
  use moose_interp_curve
  use moose_hypersurface
  class(scalar_mfunc), intent(in   ) :: Phi
  real(real64),        intent(in   ) :: x0(2), dxmax
  integer,             intent(in   ) :: dir
  integer,             intent(  out) :: istat
  type(rlist2),        intent(in   ), optional :: Xpoints
  type(hypersurf2d),   intent(in   ), optional :: boundary
  logical,             intent(in   ), optional :: reverse
  real(real64),        intent(in   ), optional :: epsabs, dxpmax, alpha
  integer,             intent(in   ), optional :: nmax
  integer,             intent(  out), optional :: iconnect
  real(real64),        intent(  out), optional :: uconnect
  type(interp_curve)                 :: this

  type(contour) :: C


  C = contour(Phi, x0, dir, dxmax, istat, Xpoints, boundary, reverse, epsabs, dxpmax, nmax, alpha, iconnect, uconnect)
  this = C%interp()

  end function interp_contour
  !-----------------------------------------------------------------------------

end module moose_contours
