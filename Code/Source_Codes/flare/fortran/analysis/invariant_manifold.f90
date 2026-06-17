module flare_invariant_manifold
  use iso_fortran_env
  implicit none


  real(real64) :: offset = 1.d-3, dx = 1.d-3

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine invariant_manifold(ix, idir, nsym, nfp, nphi, phi0, n, x, b, u)
  use moose_mpi
  use moose_math,    only: pi2
  use moose_rlist
  use moose_utils,   only: progress_bar, wait_for_all_procs
  use flare_control, only: report
  use flare_model, only: assert_equi2d, equi2d
  use flare_fieldline
  integer,      intent(in   ) :: ix, idir, nsym, nfp, nphi
  real(real64), intent(in   ) :: phi0
  integer,      intent(  out) :: n(0:nphi-1), b(0:nphi-1)
  real(real64), intent(  out) :: x(2, 0:nphi-1, 0:nfp), u(2, 0:nphi-1)

  type(fdriver) :: F
  type(rlist2)  :: tmp
  real(real64)  :: dx2, v(2,-1:1), x0(2), x1(2), y(3)
  integer       :: i, i1, istat, itrace, j, j1


  call assert_equi2d("invariant_manifold")
  itrace = -equi2d%Bt_sign * idir
  F  = fdriver()
  v  = equi2d%xpoint_stability(ix)
  x0 = equi2d%xpoint(ix)
  x1 = x0 - offset * v(:,idir)


  x = 0.d0
  u = 0.d0
  b = 0
  n = 0
  ! construct initial points at phi0
  do i=rank,nphi-1,nproc
     call F%reset()

     y(1:2) = x1
     y(3)   = phi0 - itrace * pi2/nsym * i / nphi
     istat  = F%evolve3(y, phi0)
     if (istat /= 0) call FIELDLINE_ERROR(F, istat)
     x(:,i,0) = y(1:2)
  enddo


  ! generate points on invariant manifold
  if (report) call progress_bar(0, nphi)
  do i=rank,nphi-1,nproc
     call F%reset()

     y(1:2) = x(:,i,0)
     y(3)   = phi0
     n(i)   = nfp
     do j=1,nfp
        istat = F%evolve3(y, phi0 + j*itrace*pi2/nsym)
        x(:,i,j) = y(1:2)
        if (istat == INTERSECT_BOUNDARY  .or.  istat == EDOM) then
           n(i) = j-1
           b(i) = F%nb
           u(1,i) = y(3)
           u(2,i) = F%ub(2)
           exit
        elseif (istat > 0) then
           call FIELDLINE_ERROR(F, istat)
        endif
     enddo
     if (report) call progress_bar(i+1, nphi)
  enddo
  call wait_for_all_procs(report)
  call moose_mpi_sum(x)
  call moose_mpi_sum(u)
  call moose_mpi_sum(n)
  call moose_mpi_sum(b)

  end subroutine invariant_manifold
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function footprint_size(nfp, nphi, n, x) result(w)
  use flare_model, only: equi2d
  integer,      intent(in) :: nfp, nphi, n(0:nphi-1)
  real(real64), intent(in) :: x(2, 0:nphi-1, 0:nfp)
  real(real64)             :: w

  integer :: i, j


  w = 0.d0
  do i=0,nphi-1
     if (n(i) == nfp) cycle

     j = n(i) + 1
     w = max(w, equi2d%psiN(x(:,i,j))-1.d0)
  enddo

  end function footprint_size
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine transform_u(nphi, b, u, isp, usp)
  use moose_error,    only: ERROR
  use moose_utils,    only: str
  use flare_model,    only: boundary2d
  integer,      intent(in   ) :: nphi, b(0:nphi-1), isp
  real(real64), intent(inout) :: u(2, 0:nphi-1)
  real(real64), intent(in   ) :: usp

  real(real64), allocatable :: s(:)
  real(real64) :: s0, si
  integer :: i, k


  allocate (s(0:boundary2d%P(isp)%segments()), source=boundary2d%P(isp)%accumulated_lengths())
  k = int(usp)
  s0 = s(k) + (s(k+1) - s(k)) * (usp - k)
  do i=0,nphi-1
     ! nothing to be done if there is no strike point data
     if (b(i) == 0) cycle

     ! sanity check for strike point
     if (b(i) /= isp) then
        print *, "i    = ", i
        print *, "b(i) = ", b(i)
        print *, "isp  = ", isp
        call ERROR("inconsistent boundary id ", "transform_u")
     endif

     ! convert u(2,i) into distance [cm] from equilibrium strike point usp
     k = int(u(2,i))
     u(2,i) = abs(s(k) + (s(k+1) - s(k)) * (u(2,i) - k) - s0) * 1.d2
  enddo

  end subroutine transform_u
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function sparse_output(x0, nfp, nphi, n, x, b, dmin, dmax) result(output)
  use moose_rlist
  real(real64), intent(in) :: x0(2), x(2, 0:nphi-1, 0:nfp), dmin, dmax
  integer,      intent(in) :: nfp, nphi, n(0:nphi-1), b(0:nphi-1)
  type(rlist2)             :: output

  real(real64) :: d2, dmin2, dmax2, x1(2)
  integer :: i, j


  output = rlist2()
  call output%append(x0)


  dmin2 = dmin**2
  dmax2 = dmax**2
  x1 = x0
  jloop: do j=0,nfp
  iloop: do i=0,nphi-1
     ! skip undefined points after intersection with boundary
     if (j > n(i)) cycle

     d2 = sum((x(:,i,j)-x1)**2)
     ! stop output if distance between points is too large
     if (d2 > dmax2) then
        exit jloop
     endif

     ! ignore point if distance to last point is too small
     if (d2 < dmin2) cycle

     ! save point
     x1 = x(:,i,j)
     call output%append(x1)
  enddo iloop
  enddo jloop

  end function sparse_output
  !-----------------------------------------------------------------------------

end module flare_invariant_manifold
