module flare_mmesh_base_generator
  use moose_grids,    only: qmesh
  use moose_geometry, only: curve
  use iso_fortran_env
  implicit none


  type guiding_contour_container
     class(curve), pointer :: curve
  end type guiding_contour_container


  character(len=32) :: &
     inner_boundary_method = "default"

  real(real64) :: &
     updown_symmetry_tolerance = 1.d-3, &
     min_width = -1.d0, &
     max_squeeze = 0.8d0


  logical, private :: equi2d_base_mesh = .true.

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_base_mesh_generator()
  use flare_model,            only: assert_equi2d
  use flare_mmesh_parameters, only: layout, LAYOUT_GENERIC
  use flare_mmesh_inner_boundary
  use base_mesh


  equi2d_base_mesh = layout /= LAYOUT_GENERIC
  if (equi2d_base_mesh) call assert_equi2d("init_base_mesh_generator")

  call load_inner_boundary(.true., equi2d_base_mesh)
  if (equi2d_base_mesh) call init_equi2d_base_mesh_generator()

  end subroutine init_base_mesh_generator
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function base_filename(iz) result(filename)
  use moose_utils, only: str
  integer, intent(in) :: iz
  character(:), allocatable :: filename


  filename = "base"//str(iz)//".dat"

  end function base_filename
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine generate_generic_base_mesh(iblock)
  use moose_error
  use moose_utils, only: nsubstrings, str, substring
  use moose_quantiles
  use moose_geometry,         only: loadtxt_curve
  use flare_mmesh_parameters, only: guiding_contour, tblock, radial_spacing, poloidal_spacing, &
                                    blocks, updown_symmetry, qmesh_generator
  use flare_mmesh_inner_boundary
  integer, intent(in) :: iblock

  class(qfunc), allocatable :: Qr, Qp
  type(guiding_contour_container), allocatable :: C(:)
  type(qmesh) :: this
  character(len=256) :: cmd, gc
  real(real64), allocatable :: t(:)
  integer, allocatable :: ir(:)
  integer :: i, i0, i1(-1:1), icontour, idir, ilayer, ios, n, nr, np


  ! poloidal spacing of mesh nodes
  np = tblock(iblock)%np(0)
  allocate (Qp, source=generate_quantile_function(poloidal_spacing(iblock)))
  allocate (t(0:np), source=0.d0)


  ! guiding contours & sub-layers
  gc = guiding_contour(iblock);   if (gc == "") gc = guiding_contour(-1)
  n  = nsubstrings(gc, ',')
  if (n == 0) call ERROR("missing definition of guiding contour")
  allocate (ir(-1:n), source=0)
  allocate (C(-1:n))
  C(-1)%curve => inner_boundary(0,iblock)%curve
  C(0)%curve => inner_boundary(1,iblock)%curve
  do i=1,n
     ir(i) = ir(i-1) + tblock(iblock)%nr(i-1)
     allocate (C(i)%curve, source=loadtxt_curve(substring(gc, i, ',')))
     if (.not.C(i)%curve%is_closed) call ERROR("guiding contour "//str(i)//" is not closed")
  enddo
  ir(0) = 1 ! 2nd inner boundary


  ! initialize mesh
  nr = sum(tblock(iblock)%nr(0:n-1))
  this = qmesh(nr+1, np+1)
  if (qmesh_generator == "") call ERROR("undefined qmesh_generator")
  read (qmesh_generator, *) cmd
  ! reference surface
  select case(cmd)
  case("QUASI_ORTHOGONAL")
     i0 = 1

  case("DISTANCE_GRADIENTS")
     i0 = n

  case default
     i0 = 0
  end select
  t  = C(i0)%curve%arclength_quantiles(np, Qp)
  this%x(ir(i0),:,:) = transpose(C(i0)%curve%eval(t))


  ! construct mesh layers
  i1 = [0, 0, n-1]
  if (inner_boundary_method /= "default") i1(-1) = 1
  do idir=1,-1,-2
  do i=i0,i1(idir),idir
     ilayer = min(i,i+idir)
     icontour = i+idir

     ! radial node range and spacings
     allocate (Qr, source=generate_quantile_function(radial_spacing(ilayer)))

     ! ir(i) -> ir(i+idir)
     call aux_submesh(ir(i), ir(i+idir))

     ! cleanup & prepare next layer
     call Qr%free();   deallocate (Qr)
  enddo
  enddo


  ! construct nodes on 1st inner boundary
  select case(inner_boundary_method)
  case ("default")
     ! nothing to be done here

  case ("radial_projection")
     print *, "mapping nodes to 1st boundary ..."
     call map_nodes_to_1st_boundary(this, iblock)

  case ("footpoints")
     print *, "constructing footpoints on 1st boundary ..."
     call find_footpoints_on_1st_boundary(this, iblock)

  case default
     call ERROR("invalid choice '"//trim(inner_boundary_method)//"' for inner_boundary")
  end select


  ! apply up/down symmetry
  if ((updown_symmetry < 0  .and.  iblock == 0)  .or.  &
      (updown_symmetry > 0  .and.  iblock == blocks-1)) then
     call apply_updown_symmetry(this, iblock)
  endif
  call this%savetxt(base_filename(iblock))


  ! cleanup
  call Qp%free();   deallocate (Qp)

  contains
  !.............................................................................
  subroutine aux_submesh(ir1, ir2)
  use moose_error
  use moose_utils,            only: str
  use moose_qmesh_generator
  use flare_control,          only: verbose
  integer, intent(in) :: ir1, ir2

  character(len=256) :: dummy
  type(qmesh)  :: tmp
  real(real64) :: ta, tb, ca, cb
  integer      :: m


  select case (cmd)
  ! quasi-orthogonal mesh from dummy potential between ir1 and ir2
  case ("QUASI_ORTHOGONAL")
     print *, "constructing quasi-orthogonal mesh between ", ir1, " and ", ir2, " ..."
     debug_dummyU = "dummyU_submesh"//str(ir1)//"-"//str(ir2)
     call construct_submesh(this, ir1, ir2, C(icontour)%curve, Qr, 16, min_width, max_squeeze)


  ! mesh from distance contours
  case ("DISTANCE_CONTOURS")
     print *, "constructing mesh between ", ir1, " and ", ir2, " from distance contours ..."
     call aux_qmesh_distance_contour_generator(this, ir1, ir2, C(icontour)%curve)


  ! mesh from distance gradient trace
  case ("DISTANCE_GRADIENTS")
     print *, "constructing mesh between ", ir1, " and ", ir2, " from distance gradient trace ..."
     call aux_qmesh_distance_gradient_generator(this, ir1, ir2, C(icontour)%curve)


  ! intepolate mesh in blocks between normal rays from ir1 to ir2
  case ("NORMAL_RAYS")
     print *, "constructing mesh between ", ir1, " and ", ir2, " from interpolation of rays ..."
     read  (qmesh_generator, *, iostat=ios) dummy, m;   call assert_ios()
     call normal_ray_interpolate(this, ir1, ir2, C(icontour)%curve, m, idir, Qr)


  ! interpolate between ir1 and ir2
  case ("INTERPOLATE")
     print *, "interpolating mesh between ", ir1, " and ", ir2, " ..."
     !this%x(ir2,:,:) = transpose(C(icontour)%curve%arclength_discretization(np+1, Qp))

     ! update t for curve domain
     ta = t(lbound(t,1));   ca = C(icontour)%curve%a
     tb = t(ubound(t,1));   cb = C(icontour)%curve%b
     t = ca  +  (cb-ca) * (t-ta) / (tb-ta)
     t = min(max(t, ca), cb)   ! truncate values (finite accuracy!)

     this%x(ir2,:,:) = transpose(C(icontour)%curve%eval(t))
     call interpolate_submesh(this, ir1, ir2, 1, Qr)


  ! invalid generator
  case default
     call ERROR("invalid qmesh generator '"//trim(cmd)//"'")
  end select


  ! save submesh
  if (verbose) then
     tmp = qmesh(this%u(min(ir1,ir2):max(ir1,ir2),:), this%v(min(ir1,ir2):max(ir1,ir2),:))
     call tmp%savetxt("submesh"//str(ir1)//"-"//str(ir2))
  endif

  end subroutine aux_submesh
  !.............................................................................
  subroutine assert_ios()


  if (ios == 0) return
  call ERROR("missing or invalid parameter for "//trim(cmd))

  end subroutine assert_ios
  !.............................................................................
  end subroutine generate_generic_base_mesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine apply_updown_symmetry(this, iblock)
  use moose_mpi
  use moose_error
  use moose_math,             only: pi
  use flare_model,            only: bfield
  use flare_mmesh_parameters, only: T
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: iblock

  real(real64) :: d, dmax, phi, x(2), x0(2)
  integer :: ir, ip, jp, nr, np


  phi  = T(iblock)%phi_base / 180.d0 * pi
  x0   = bfield%equi%magnetic_axis(phi)
  dmax = 0.d0
  nr   = this%n(1)-1
  np   = this%n(2)-1
  do ir=0,nr
  do ip=0,np/2
     jp = np-ip
     x(1) = (this%x(ir,ip,1) + this%x(ir,jp,1)) / 2
     x(2) = (this%x(ir,ip,2) - this%x(ir,jp,2)) / 2
     d    = sqrt(sum((this%x(ir,ip,:) - x)**2)) / x0(1)
     if (d > updown_symmetry_tolerance) then
        print *, "ir, ip, jp = ", ir, ip, jp
        print *, "x(ir,ip)   = ", this%x(ir,ip,:)
        print *, "x(ir,jp)   = ", this%x(ir,jp,:)
        print *, "d          = ", d
        call this%savetxt("ERROR_UPDOWN_SYMMETRY")
        call ERROR("up/down symmetry exceeds tolerance")
     endif
     dmax = max(d,dmax)
     this%x(ir,ip,:) = x
     this%x(ir,jp,:) = [x(1), -x(2)]
  enddo
  enddo
  if (rank == 0) print *, "up/down symmetry deviation: ", dmax

  end subroutine apply_updown_symmetry
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine generate_base_mesh(iblock)
  use base_mesh
  integer, intent(in) :: iblock

  type(qmesh) :: M
  logical :: usf
  integer :: nr, np


  print 1000, iblock
 1000 format(3x,"- Generating base mesh for block ",i0)


  if (equi2d_base_mesh) then
     call generate_equi2d_base_mesh(iblock)
  else
     call generate_generic_base_mesh(iblock)
  endif

  end subroutine generate_base_mesh
  !-----------------------------------------------------------------------------

end module flare_mmesh_base_generator
