module flare_mmesh_inner_boundary
  use iso_fortran_env
  use moose_geometry, only: curve, hypersurf2d
  use flare_tasks
  implicit none


  ! number of linear segments for approximation
  integer :: &
     hypersurf2d_segments = 1024



  type boundary_data
     ! actual geometry
     class(curve), pointer :: curve

     ! hypersurf2d representation
     type(hypersurf2d), pointer :: hypersurf2d

     ! mean, variation, minimum and maximum of normalized poloidal flux
     real(real64) :: psiN, delta_psiN, min_psiN, max_psiN
  end type boundary_data
  type(boundary_data), allocatable :: inner_boundary(:,:)

  real(real64) :: psiN_1

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function inner_boundary_filename(ir, iblock, suffix) result(filename)
  use moose_utils, only: ordinal, str
  integer,          intent(in) :: ir, iblock
  character(len=*), intent(in) :: suffix
  character(:), allocatable    :: filename


  filename = ordinal(ir+1) // "_inner_boundary" // str(iblock) // "." // suffix

  end function inner_boundary_filename
  !-----------------------------------------------------------------------------
  function points_filename(ir, iblock) result(filename)
  use moose_utils, only: ordinal
  integer,       intent(in) :: ir, iblock
  character(:), allocatable :: filename


  filename = inner_boundary_filename(ir, iblock, "txt")

  end function points_filename
  !-----------------------------------------------------------------------------
  function curve_filename(ir, iblock) result(filename)
  use moose_utils, only: ordinal
  integer,       intent(in) :: ir, iblock
  character(:), allocatable :: filename


  filename = inner_boundary_filename(ir, iblock, "dat")

  end function curve_filename
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine make_points_default(p1, p2, nsymmetry, phi_start, nblocks, npoints)
  !
  ! Construct Poincare maps for inner simulation boundary (default layout).
  !
  ! **Parameters:**
  !
  ! :p1, p2:   Reference points on first and second flux surfaces at inner simulation boundary in cylindrical coordinates.
  !
  use moose_mpi
  use moose_math, only: pi2
  use flare_control
  use flare_poincare_map
  use flare_mmesh, only: sector
  real(real64), intent(in) :: p1(3), p2(3), phi_start
  integer,      intent(in) :: nsymmetry, nblocks, npoints

  type(poincare_map), allocatable :: P(:,:)
  real(real64) :: dphi, p0(3,0:1), phiX
  integer      :: i, j, n


  ! print task parameters
  call begin_task()
  if (rank == 0) then
     print *, "Constructing Poincare maps for inner simulation boundary ..."
     print *
     print 1001
     print 1002, 1, p1(1:2), p1(3)
     print 1002, 2, p2(1:2), p2(3)
     print *

     dphi = sector(nsymmetry)
     phiX = phi_start + dphi/2/nblocks
     print 1010
     do j=0,nblocks-1
        print 1011, j, phiX + j * dphi/nblocks
     enddo
     phiX = phiX / 360.d0 * pi2
  endif
  call proc(0)%broadcast(phiX)
  p0(:,0) = p1;   p0(:,1) = p2;   p0(3,:) = p0(3,:) / 360.d0 * pi2
 1001 format(3x,"- Reference points (R[m], Z[m], phi[deg]) for flux surfaces")
 1002 format(8x,i0,":",3x,f0.3,", ",f0.3,", ",f0.3)
 1010 format(3x,"- Default layout:  base grid location [deg]")
 1011 format(16x,i4,16x,f8.3)


  ! construct Poincare maps
  n = nblocks;   if (nsymmetry < 0) n = 2*nblocks
  allocate (P(0:n-1, 0:1))
  do i=rank,1,nproc
     P(:,i) = poincare_maps(p0(:,i), 1, phiX, abs(nsymmetry), npoints, n)

     do j=0,nblocks-1
        call P(j,i)%savetxt(points_filename(i, j))
     enddo
  enddo


  ! cleanup
  call finalize_task()
  do i=rank,1,nproc
     do j=0,n-1
        call P(j,i)%free()
     enddo
  enddo
  deallocate (P)

  end subroutine make_points_default
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine make_points_usr(p1, p2, symmetry, phi_base, npoints)
  !
  ! Construct Poincare maps for inner simulation boundary (user defined base positions).
  !
  ! **Parameters:**
  !
  ! :p1, p2:   Reference points on first and second flux surfaces at inner simulation boundary in cylindrical coordinates.
  !
  use moose_mpi
  use moose_math, only: pi2
  use flare_control
  use flare_poincare_map
  real(real64), intent(in) :: p1(3), p2(3), phi_base(:)
  integer,      intent(in) :: symmetry, npoints

  type(poincare_map), allocatable :: P(:,:,:)
  real(real64) :: p0(3,0:1), phiX
  integer      :: i, j, k, nblocks


  ! print task parameters
  call begin_task()
  nblocks = size(phi_base)
  if (rank == 0) then
     print *, "Constructing Poincare maps for inner simulation boundary ..."
     print *
     print 1001
     print 1002, 1, p1(1:2), p1(3)
     print 1002, 2, p2(1:2), p2(3)

     print *
     print 1010
     do j=0,nblocks-1
        print 1011, j, phi_base(j+1)
     enddo
  endif
  p0(:,0) = p1;   p0(:,1) = p2;   p0(3,:) = p0(3,:) / 360.d0 * pi2
 1001 format(3x,"- Reference points (R[m], Z[m], phi[deg]) for flux surfaces")
 1002 format(8x,i0,":",3x,f0.3,", ",f0.3,", ",f0.3)
 1010 format(3x,"- Non-default layout:  base grid location [deg]")
 1011 format(20x,i4,16x,f8.3)


  ! construct Poincare maps
  allocate (P(1,0:nblocks-1, 0:1))
  k = -1
  do i=0,1
     do j=0,nblocks-1
        k = k + 1
        if (.not.mod(k,nproc) == rank) cycle

        phiX = phi_base(j+1) / 360.d0 * pi2
        P(:,j,i) = poincare_maps(p0(:,i), 1, phiX, symmetry, npoints, 1)
        call P(1,j,i)%savetxt(points_filename(i, j))
        call P(1,j,i)%free()
     enddo
  enddo


  ! cleanup
  call finalize_task()
  deallocate (P)

  end subroutine make_points_usr
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bspline_multifit(nblocks, spline_order, nctrl, epsabs, fit_method, lambda1, lambda2)
  !
  ! Fit B-Spline curve to Poincare maps for inner simulation boundary.
  !
  use moose_mpi
  use moose_math,     only: pi
  use moose_bfit
  use flare_model
  use flare_equi2d, only: equi2d_bfield => equi2d
  use flare_poincare_map
  use flare_control
  integer,      intent(in) :: nblocks, spline_order, nctrl, fit_method
  real(real64), intent(in) :: epsabs, lambda1, lambda2

  type(poincare_map) :: P
  type(bfit) :: B
  real(real64) :: theta, chisq_dof
  logical :: is_equi2d
  integer :: i, j, k


  is_equi2d = .false.
  select type(equi => bfield%equi)
  class is(equi2d_bfield)
     is_equi2d = .true.
  end select


  ! print task parameters
  call begin_task()
  if (rank == 0) then
     print *, "Constructing B-Spline fit to Poincare maps for inner simulation boundary ..."
     print *
     print *, "    Surface     Toroidal position [deg]      fit coeffs.    chisq/dof"
  endif


  k = 0
  do i=0,1
     do j=0,nblocks-1
        k = 2*j + i;   if (mod(k,nproc) /= rank) cycle

        ! load poincare map
        P = poincare_map(points_filename(i, j))

        ! fit B-Spline
        if (nctrl == 0) then
           B = P%bspline_autofit(spline_order, epsabs, knot_balancing=.true., method=fit_method)
        else
           B = P%bspline_multifit(nctrl, spline_order, 0.d0, .true., fit_method, epsabs, lambda1, lambda2)
        endif
        chisq_dof = B%e(0) / B%nctrl
        call B%savetxt(curve_filename(i, j))
        print 1000, i, 180*P%phiX/pi, B%nctrl, chisq_dof

        ! output psiN(theta) contour
        if (is_equi2d) call theta_psiN_contour(inner_boundary_filename(i, j, "plt"))
     enddo
  enddo
 1000 format(10x,i1,12x,f10.3,12x,i8,8x,g14.7)


  call mpi_barrier_world()
  call finalize_task()

  contains
  !.............................................................................
  subroutine  theta_psiN_contour(filename)
  use moose_txtio, only: savetxt
  character(len=*), intent(in) :: filename

  real(real64), allocatable :: x(:,:)
  real(real64) :: p(2), theta, psiN
  integer :: i, n


  n = hypersurf2d_segments
  allocate (x(2, 0:n))
  do i=0,n
     p = B%eval(i, n)
     theta = bfield%equi%poloidal_angle(p) / pi * 180
     psiN = bfield%equi%psiN(p)
     x(1, i) = theta
     x(2, i) = psiN
  enddo

  do i=1,n
     if (x(1, i) - x(1, i-1) < 0.d0) x(1, i:) = x(1, i:) + 360.d0
  enddo
  call savetxt(filename, transpose(x))

  end subroutine  theta_psiN_contour
  !.............................................................................
  end subroutine bspline_multifit
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load_inner_boundary(hypersurf2d, equi2d_psiN)
  use moose_error
  use moose_utils,            only: ordinal
  use moose_geometry,         only: loadtxt_curve, polygon2d, hypersurf2d_approximation
  use flare_mmesh_parameters, only: blocks, nr_perturbed
  logical, intent(in) :: hypersurf2d, equi2d_psiN

  integer :: i, j


  print 1000
  if (equi2d_psiN) then
     print 1001
  else
     print 1002
  endif
 1000 format(3x,'- Inner simulation boundaries:')
 1001 format(8x,'block #:   mean(psiN),              min(psiN),  max(psiN)')
 1002 format(8x,'block #')


  psiN_1 = 0.d0
  allocate (inner_boundary(0:nr_perturbed-1, 0:blocks-1))
  do i=0,blocks-1
     write (6, 1010, advance='no') i
     do j=0,nr_perturbed-1
        inner_boundary(j,i) = load_boundary_data(curve_filename(j,i))
     enddo
     psiN_1 = psiN_1 + inner_boundary(nr_perturbed-1,i)%psiN
     ! check for intersections between boundary contours (due to inappropriate representation)
     if (nr_perturbed > 1) then
        associate (S => inner_boundary(nr_perturbed-1,i)%hypersurf2d%P(1))
        if (S%intersects(inner_boundary(0,i)%hypersurf2d%P(1))) then
           print *
           call ERROR("intersection between inner boundary contours detected")
        endif
        end associate
     endif

     if (equi2d_psiN) then
        associate (S => inner_boundary(nr_perturbed-1,i))
        write (6, 1011) S%psiN, S%delta_psiN, S%min_psiN, S%max_psiN

        ! assert that inner boundary is inside separatrix (psiN = 1)
        if (S%max_psiN > 1.d0) call ERROR("inner boundary appears to be outside of separatrix")
        end associate
     else
        write (6, 1012)
     endif
  enddo
  psiN_1 = psiN_1 / blocks
  print *
 1010 format(8x,i7)
 1011 format(4x,f8.6,' +/- ',f8.6,2(4x,f8.6))
 1012 format()

  contains
  !.............................................................................
  function load_boundary_data(filename) result(this)
  use flare_model, only: equi2d
  character(len=*), intent(in) :: filename
  type(boundary_data)          :: this

  real(real64) :: psiN, delta
  integer :: i, n


  allocate (this%curve, source=loadtxt_curve(filename))


  ! set up hypersurf2d approximation of curve geometry
  if (hypersurf2d  .or.  equi2d_psiN) then
     allocate (this%hypersurf2d, source=hypersurf2d_approximation(this%curve, hypersurf2d_segments))
  endif


  this%psiN       = 0.d0
  this%delta_psiN = 0.d0
  this%min_psiN   = 1.d0
  this%max_psiN   = 0.d0
  ! set up normalized poloidal flux parameters
  if (equi2d_psiN) then
     n = this%hypersurf2d%P(1)%segments()
     if (n == 0) call ERROR("discretization of inner boundary failed")
     do i=1,n
        psiN            = equi2d%psiN(this%hypersurf2d%P(1)%node(i))
        delta           = psiN - this%psiN
        this%psiN       = this%psiN + delta/i
        this%delta_psiN = this%delta_psiN + delta*(psiN - this%psiN)
        this%min_psiN   = min(this%min_psiN, psiN)
        this%max_psiN   = max(this%max_psiN, psiN)
     enddo
     this%delta_psiN = sqrt(this%delta_psiN / (n - 1.d0))
  endif

  end function load_boundary_data
  !.............................................................................
  end subroutine load_inner_boundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine map_nodes_to_1st_boundary(this, iblock)
  !
  ! set up nodes on 1st inner boundary by extrapolating from the 2nd inner boundary
  ! along the radial direction of the mesh (given by the next node in radial direction)
  !
  use moose_error
  use moose_utils, only: str
  use moose_grids, only: qmesh
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: iblock

  type(qmesh)  :: M
  real(real64) :: x1(2), x2(2), x0(2), x(2), t, u(1)
  integer :: iu, j, n


  do j=0,this%n(2)-1
     x1 = this%x(1,j,:)
     x2 = this%x(2,j,:)
     x0 = x1 + 10*(x1-x2)

     if (inner_boundary(0,iblock)%hypersurf2d%intersect(x1, x0, x, t, n, u)) then
        this%x(0,j,:) = x
     else
        ! save direction for node extrapolation
        open  (newunit=iu, file="NODE_EXTRAPOLATION")
        write (iu, *) x2
        write (iu, *) x1
        write (iu, *) x0
        close (iu)

        ! save mesh at 2nd inner boundary
        M = qmesh(this%x(1:2,:,1), this%x(1:2,:,2))
        call M%savetxt("MESH_AT_2ND_INNER_BOUNDARY"//str(iblock))

        print *
        print 9000
        print 9001, iblock, iblock
        print *
        call ERROR("mapping from 2nd to 1st boundary failed at poloidal index j = "//str(j))
     endif
  enddo
 9000 format("run the following command to verify the given geometry:"/)
 9001 format("mview MESH_AT_2ND_INNER_BOUNDARY",i0," -mplot NODE_EXTRAPOLATION,ko- 1st_inner_boundary",i0,".dat,r")

  end subroutine map_nodes_to_1st_boundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine find_footpoints_on_1st_boundary(this, iblock)
  !
  ! set up nodes on 1st inner boundary as footpoints from mesh nodes along 2nd inner boundary
  !
  use moose_error
  use moose_utils, only: str
  use moose_grids, only: qmesh
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: iblock

  real(real64), allocatable :: P(:,:), t(:), e(:)
  integer :: ierr, n


  n = size(this%x, 2)
  allocate (P(2,0:n-1), t(0:n-1), e(0:n-1))
  call inner_boundary(0,iblock)%curve%find_footpoints(n, transpose(this%x(1,:,:)), t, P, e, ierr)
  if (ierr /= 0) call ERROR("footpoint construction on inner boundary failed")
  this%x(0,:,:) = transpose(P)

  end subroutine find_footpoints_on_1st_boundary
  !-----------------------------------------------------------------------------

end module flare_mmesh_inner_boundary
