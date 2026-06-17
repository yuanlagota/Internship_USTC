subroutine flare_task_invariant_manifold(ix, idir, nsym, nfp, nphi, phi0, dmin, dmax, output)
  !
  ! Construct invariant manifold (perturbed separatrix) for selected X-point.
  !
  ! **Parameters:**
  !
  ! :ix:       X-point number.
  !
  ! :idir:     Direction of manifold (1: forward, -1: backward).
  !
  ! :nsym:     Toroidal symmetry number.
  !
  ! :nfp:      Max. number for field periods for tracing.
  !
  ! :nphi:     Sample resolution.
  !
  ! :phi0:     Reference position along toroidal direction [deg].
  !
  ! :dmin:     Min. distance between points on contour for sparse output (0: raw output).
  !
  ! :dmax:     Max. distance between points on contour (truncate contour for smooth output).
  !
  ! :output:   Filename for output of invariant manifold contour.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math,  only: pi
  use moose_rlist
  use moose_contours, only: xcontour, XCONTOUR_BRANCH, XCONTOUR_POSITIVE_ORIENTATION
  use flare_model, only: assert_equi2d, equi2d
  use flare_fluxsurf2d, only: separatrix2d
  use flare_invariant_manifold
  use flare_tasks
  implicit none
  integer,          intent(in) :: ix, idir, nsym, nfp, nphi
  real(real64),     intent(in) :: phi0, dmin, dmax
  character(len=*), intent(in) :: output

  real(real64), allocatable :: x(:,:,:), u(:,:)
  integer,      allocatable :: b(:), n(:)
  type(xcontour) :: sepx
  type(rlist2) :: M
  real(real64) :: x0(2), xsp(2), s, usp
  integer      :: i, isp, iu, j, k


  call begin_task()
  call assert_equi2d("invarian_manifolds")
  x0 = equi2d%xpoint(ix)
  if (report) then
     print 1000, ix, x0
     print 1001, phi0
     print 1002, nsym

     ! find strike points of separatrix
     sepx = separatrix2d(ix)
     k = XCONTOUR_BRANCH(idir, XCONTOUR_POSITIVE_ORIENTATION)
     i = -1;   if (idir == 1) i = 0
     xsp = sepx%branch(k)%point(i)
     isp = -sepx%iconnect(k)
     usp = sepx%uconnect(k)
     print 1003, xsp, isp, usp

     print 1004, nfp, nphi
  endif
  call proc(0)%broadcast(isp)
  call proc(0)%broadcast(usp)
 1000 format(1x,"Construcing invariant manifolds for X-point ",i0,": (",f0.3,", ",f0.3,") m",/)
 1001 format(3x,"- Toroidal position for R-Z slice: ",f0.3," deg",/)
 1002 format(3x,"- Toroidal symmetry: ",i0,/)
 1003 format(3x,"- Strike point:  (",f0.3,", ",f0.3,") m on boundary ",i0," at u = ",f0.3,/)
 1004 format(3x,"- Tracing ",i0," field periods with ",i0," sample points",/)


  allocate (n(0:nphi-1), x(2, 0:nphi-1, 0:nfp), b(0:nphi-1), u(2,0:nphi-1))
  call invariant_manifold(ix, idir, nsym, nfp, nphi, phi0/180.0*pi, n, x, b, u)
  call transform_u(nphi, b, u, isp, usp)
  if (rank == 0) then
     s = footprint_size(nfp, nphi, n, x)
     print 2000, s

     if (dmin == 0.d0) then
        open  (newunit=iu, file=output)
        write (iu, *) nsym, nfp, nphi, s
        write (iu, *) n
        write (iu, *) b
        write (iu, *) x
        write (iu, *) u
        close (iu)
     else
        M = sparse_output(x0, nfp, nphi, n, x, b, dmin, dmax)
        call M%savetxt(output)
     endif
  endif
 2000 format(/,3x,"- Footprint size: ",f0.4," [psiN]")


  deallocate (n, x, b, u)
  call finalize_task()

end subroutine flare_task_invariant_manifold
