module flare_boundary
  use iso_fortran_env
  implicit none


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
!  function firstwall(phi)
!  !
!  ! construct firstwall torosurf with cross-sections sampled at phi
! NOTE: first wall may be split into 2 stellarator symmetric parts
!  !
!  use moose_geometry, only: polygon2d, torosurf
!  use moose_torosurf, only: setup_torosurf
!  use flare_model,    only: boundary
!  real(real64), intent(in) :: phi(:)
!  type(torosurf)           :: firstwall
!
!  type(polygon2d) :: rzslice
!  logical :: axisymmetric_firstwall
!  integer :: iboundary, iphi, nsym


  ! determined boundary index and symmetry of firstwall
!  call aux_firstwall_rzslice(phi(1), iboundary, rzslice)
!  axisymmetric_firstwall = iboundary <= boundary%nA
!  if (axisymmetric_firstwall) then
!     nsym = boundary%A(iboundary)%symmetry
!  else
!     nsym = boundary%T(iboundary - boundary%nA)%nsym
!  endif


!  ! generate contours of firstwall
!  firstwall = torosurf(size(phi)-1, rzslice%segments(), nsym)
!  firstwall%phi = phi
!  if (axisymmetric_firstwall) then
!     firstwall%rz(:,:,:) = spread(transpose(boundary%A(iboundary)%P%nodes()), 3, size(phi))
!  else
!     do iphi=1,size(phi)
!        firstwall%rz(:,:,iphi-1) = boundary%T(iboundary - boundary%nA)%slice(phi(iphi))
!     enddo
!  endif
!  call setup_torosurf(firstwall)

!  end function firstwall
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function firstwall_rzslice(phi) result(rzslice)
  !
  ! slice through first wall at *phi* [rad]
  !
  ! first wall is:
  !   1) poloidally closed
  !   2) contains magnetic axis at phi
  !   3) is the one with the smallest cross-section at phi that matches 1) and 2)
  !
  use moose_geometry, only: polygon2d
  real(real64), intent(in) :: phi
  type(polygon2d)          :: rzslice

  integer :: iboundary


  call aux_firstwall_rzslice(phi, iboundary, rzslice)

  end function firstwall_rzslice
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_firstwall_rzslice(phi, iboundary, rzslice)
  use moose_error
  use moose_math
  use moose_geometry, only: polygon2d, axisurf, torosurf
  use flare_model,    only: boundary, bfield
  real(real64),    intent(in   ) :: phi
  integer,         intent(  out) :: iboundary
  type(polygon2d), intent(  out) :: rzslice

  real(real64) :: minA, x(2), dphi, dphi2, phiL, phiU
  integer :: i


  x = bfield%equi%magnetic_axis(phi)
  minA = huge(1.d0)

  do i=1,boundary%nsurfaces
  select type(S => boundary%surfaces(i)%geometry)
  ! axisymmetric surfaces
  type is (axisurf)
     call check(i, S%P)

  ! non-axisymmetric surfaces
  type is (torosurf)
     if (.not.S%includes(phi, .true.)) cycle

     phiL = S%phi(0)
     phiU = S%phi(S%nu)
     dphi = pi2 / S%symmetry
     dphi2 = dphi / 2
     if (phi < dphi2) then
        ! lower boundary must be at 0
        if (phiL /= 0.d0) cycle

        ! upper boundary must be at dphi or dphi/2
        if (phiU /= dphi  .and.  phiU /= dphi2) cycle

     else
        ! lower boundary must be at 0 or dphi (stellarator symmetry)
        if (phiL /= 0.d0  .and.  phiL /= dphi) cycle

        ! upper boundary must be at dphi or phi2 (stellarator symmetry)
        if (phiU /= dphi  .and.  phiU /= dphi2) cycle
     endif

     call check(i, S%polygon2d(phi))
  end select
  enddo


  if (minA == huge(1.d0)) call ERROR("cannot determine first wall")


  contains
  !.............................................................................
  subroutine check(i, P)
  integer,         intent(in) :: i
  type(polygon2d), intent(in) :: P

  real(real64) :: A


  if (P%is_closed()  .and.  P%winding_number(x) /= 0) then
     A = P%area()
     if (A < minA) then
        minA = A
        iboundary = i
        rzslice = P
     endif
  endif

  end subroutine check
  !.............................................................................
  end subroutine aux_firstwall_rzslice
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine convert_coordinates(i, u1, u2)
  !
  ! convert mesh coordinates (u1, u2) to surface coordinates (phi, v) for i-th boundary
  !
  use flare_model, only: assert_model, boundary
  integer,      intent(in   ) :: i
  real(real64), intent(inout) :: u1(:), u2(size(u1))

  real(real64) :: vphi(2)
  integer :: k


  call assert_model()
  do k=1,size(u1)
     vphi = boundary%vphi(i, [u1(k), u2(k)])
     u1(k) = vphi(2)
     u2(k) = vphi(1)
  enddo

  end subroutine convert_coordinates
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function strike_point_density(F, dphi, ds, report) result(this)
  !
  ! Strike points on PFCs (from *fieldline_connection* task) are gathered in a
  ! surface mesh with resolution *dphi* [deg] and *ds* [m].
  !
  use moose_error
  use moose_utils
  use moose_math, only: pi
  use moose_units
  use moose_r3grid
  use moose_geometry, only: axisurf, torosurf
  use moose_hypermesh3d
  use moose_dataset
  use flare_version
  use flare_model, only: boundary
  type(dataset), intent(in) :: F
  real(real64),  intent(in) :: dphi, ds
  logical,       intent(in), optional :: report
  type(dataset)             :: this

  type(hypermesh3d), pointer :: M
  real(real64), pointer :: x(:)
  real(real64) :: u(2), f100, maxp
  integer, allocatable :: summary(:)
  integer :: i, k, k1, k2, nbwd, nfwd, n, nloss, npoints


  nbwd = F%data_index("n_bwd")
  nfwd = F%data_index("n_fwd")
  if (nbwd == 0  .or.  nfwd == 0) then
     call ERROR("dataset F does not include strike point data")
  endif


  ! initialize data set
  allocate (M, source = hypermesh3d(boundary, dphi, ds))
  this = dataset(5, M, CELL_DATA, kind="strike_point_density")
  call this%set_metadata(1, 'count_bwd', "Strike point count (backward direction)")
  call this%set_metadata(2, 'count_fwd', "Strike point count (forward direction)")
  call this%set_metadata(3, 'area',      "Surface area per cell", "m**2")
  call this%set_metadata(4, 'p_bwd',     "Strike point density (backward direction)", "m**(-2)")
  call this%set_metadata(5, 'p_fwd',     "Strike point density (forward direction)", "m**(-2)")
  call this%set_expression("count", "count_bwd + count_fwd", "Strike point count")
  call this%set_expression("p",     "p_bwd + p_fwd",         "Strike point density", "m**(-2)")
  call this%annotations%set("flare_version", version)


  ! count field line strike points on boundary surfaces
  allocate (summary(size(M%refined_tpzmesh3d)), source = 0)
  do k=0,F%npoints-1
     x => F%element(k)

     ! backward strike point
     n = int(x(nbwd))
     if (n /= 0) then
        u = x(nbwd+1:nbwd+2)
        i = M%cell_index(n, u)
        this%values(1,i) = this%values(1,i) + 1.d0 / boundary%surfaces(n)%geometry%symmetry
        summary(n) = summary(n) + 1
     endif

     ! forward strike point
     n = int(x(nfwd))
     if (n /= 0) then
        u = x(nfwd+1:nfwd+2)
        i = M%cell_index(n, u)
        this%values(2,i) = this%values(2,i) + 1.d0 / boundary%surfaces(n)%geometry%symmetry
        summary(n) = summary(n) + 1
     endif
  enddo


  ! compute strike point density
  npoints = sum(summary)
  nloss = count(F%values(nbwd,:) == 0.d0) + count(F%values(nfwd,:) == 0.d0)
  this%values(3,:) = M%area()
  this%values(4,:) = this%values(1,:) / (2*F%npoints - nloss) / this%values(3,:)
  this%values(5,:) = this%values(2,:) / (2*F%npoints - nloss) / this%values(3,:)


  ! report results
  if (logical_option(.false., report)) then
     print 1000
     print 1001, 2*F%npoints
     print 1002, npoints
     print 1003, nloss
     print *
     print 1004
     print *, "      ----------------------------------------------------------"
     do n=1,size(M%refined_tpzmesh3d)
        k1 = M%cell_offset(n)
        k2 = M%cell_offset(n+1) - 1
        f100 = 1.d2 * summary(n) / npoints
        maxp = maxval(sum(this%values(4:5,k1:k2), dim=1))
        print 1005, boundary%surfaces(n)%key, f100, maxp
     enddo
     print *
  endif
 1000 format(3x,"- Summary:")
 1001 format(8x,"Total number of field lines (fwd + bwd):     ", i8)
 1002 format(8x,"Number of strike points:                     ", i8)
 1003 format(8x,"Number of truncated field lines:             ", i8)
 1004 format(33x,"contribution         max. value [m**(-2)]")
 1005 format(1x,a32,2x,f8.3," %",8x,e12.5)

  end function strike_point_density
  !-----------------------------------------------------------------------------

end module flare_boundary
