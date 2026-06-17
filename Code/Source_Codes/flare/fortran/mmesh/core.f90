module flare_mmesh_core
  use iso_fortran_env
  use flare_mmesh
  use flare_mmesh_parameters, only: symmetry, updown_symmetry, blocks, T
  use flare_mmesh_layout,     only: layers
  implicit none

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine construct_rpath2d(iz, psiN0)
  !
  ! Construct core domain in zone *iz* from paths along grad psiN of toroidally
  ! symmetric equilibrium.
  !
  use moose_mpi
  use moose_error
  use moose_math,  only: linspace
  use flare_control
  use flare_model, only: assert_equi2d, equi2d, bfield
  use flare_rpath2d
  integer,      intent(in) :: iz
  real(real64), intent(in) :: psiN0

  type(rpath2d_driver) :: R
  real(real64), allocatable :: xtmp(:,:,:,:), psiNgrid(:)
  real(real64) :: x(2), psiN
  integer :: i, ir, ir1, ip, ip1, ip2, it, istat, n, nt


  if (rank == 0) then
     print 1000
  endif
 1000 format(3x,"- Tracing paths along grad psiN ..."/)


  ! allocate local workspace
  ir1 = workspace(iz)%r_surf_pl_trans_range(1)
  ip1 = workspace(iz)%p_surf_pl_trans_range(1)
  ip2 = workspace(iz)%p_surf_pl_trans_range(2)
  nt  = workspace(iz)%n(3)
  allocate (xtmp(0:ir1-1, ip1:ip2, 0:nt-1, 3), source=0.d0)
  allocate (psiNgrid(0:ir1))


  ! construct paths along grad psiN
  R = rpath2d_driver(RPATH2D_PSIN)
  i = 0
  n = nt * (ip2-ip1+1)
  call progress_bar(0, n)
  do it=0,nt-1
  do ip=ip1,ip2
     i = i + 1;   if (mod(i, nproc) /= rank) cycle
     if (verbose) print *, it, ip
     call R%reset()

     ! set reference point
     x(1) = workspace(iz)%r(ir1, ip, it)
     x(2) = workspace(iz)%z(ir1, ip, it)
     psiN = equi2d%psiN(x)
     psiNgrid = linspace(psiN0, psiN, ir1+1)

     ! step towards psiN0
     do ir=ir1-1,0,-1
        istat = R%evolve(psiN, psiNgrid(ir), x)
        if (istat /= 0) call ERROR("rpath2d_driver%evolve_apply failed", "construct_rpath2d")
        xtmp(ir,ip,it,1:2) = x
        xtmp(ir,ip,it,3)   = bfield%bmod([x(1), x(2), workspace(iz)%phi(it)])
     enddo
     call progress_bar(i, n)
  enddo
  enddo
  call finalize_progress_bar()
  

  ! collect results
  call moose_mpi_sum(xtmp)
  workspace(iz)%r(0:ir1-1, ip1:ip2, 0:nt-1) = xtmp(:,:,:,1)
  workspace(iz)%z(0:ir1-1, ip1:ip2, 0:nt-1) = xtmp(:,:,:,2)
  workspace(iz)%b(0:ir1-1, ip1:ip2, 0:nt-1) = xtmp(:,:,:,3)


  ! cleanup
  call R%free()
  deallocate (xtmp, psiNgrid)

  end subroutine construct_rpath2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine construct_flux_surfaces(iz, a)
  !
  ! Construct core domain in zone *iz* from flux surfaces.
  !
  use moose_mpi
  use moose_utils,   only: str
  use moose_math,    only: linspace, pi
  use moose_grids
  use moose_quantiles
  use flare_model
  use flare_fluxsurf3d
  use flare_control, only: verbose
  use flare_mmesh_base_generator, only: apply_updown_symmetry
  integer,      intent(in) :: iz
  real(real64), intent(in) :: a(0:workspace(iz)%r_surf_pl_trans_range(1)-1)

  type(qmesh) :: base
  type(fluxsurf3d) :: F
  type(interp_qfunc) :: Qp
  real(real64) :: dr, dz, r(3), r0(3), r1(3), theta, arcl(0:workspace(iz)%n(2)-1)
  integer :: iblock, ir, ir1, ip, it, np


  iblock = iz / layers
  it     = T(iblock)%it_base
  ir1    = workspace(iz)%r_surf_pl_trans_range(1)
  np     = workspace(iz)%n(2)
  if (rank == 0) then
     print 1000, workspace(iz)%phi(it) / pi * 180.d0
  endif
  base = qmesh(ir1, workspace(iz)%n(2))
 1000 format(3x,"- Constructing flux surface contours at ",f0.3," deg"/)


  ! magnetic axis
  r0(3)   = workspace(iz)%phi(it)
  r0(1:2) = bfield%equi%magnetic_axis(r0(3))
  ! reference point on EMC3 boundary
  r1(1) = workspace(iz)%r(ir1,0,it)
  r1(2) = workspace(iz)%z(ir1,0,it)
  r1(3) = r0(3)
  theta = atan2(r1(2)-r0(2), r1(1)-r0(1)) / pi * 180.d0
  ! poloidal spacing
  arcl = 0.d0
  do ip=1,np-1
     dr = workspace(iz)%r(ir1,ip,it) - workspace(iz)%r(ir1,ip-1,it)
     dz = workspace(iz)%z(ir1,ip,it) - workspace(iz)%z(ir1,ip-1,it)
     arcl(ip) = arcl(ip-1) + sqrt(dr**2 + dz**2)
  enddo
  Qp = interp_cdf(arcl, linspace(0.d0, 1.d0, np))


  ! construct flux surface contour & discretization based on Qp
  do ir=rank,ir1-1,nproc
     r = r0 + a(ir) * (r1-r0)
     print *, ir, r(1:2)

     F = fluxsurf3d(r, symmetry, 1, phi0=r0(3), theta0=theta, param='arclength')
     base%x(ir,:,:) = transpose(F%slice(0)%discretization(np, Qp))
  enddo
  call base%mpi_sum()
  ! apply up/down symmetry (if necessary)
  if ((updown_symmetry < 0  .and.  iblock == 0)  .or.  &
      (updown_symmetry > 0  .and.  iblock == blocks-1)) then
     call apply_updown_symmetry(base, iblock)
  endif

  if (verbose  .and.  rank == 0) call base%savetxt("core_base"//str(iblock)//".dat")


  ! trace contours
  call mpi_barrier_world()
  if (rank == 0) then
     print *
     print 2000
  endif
  call aux_trace_qmesh(iz, base)
 2000 format(3x,"- Tracing flux surface contours for 3-D mesh ..."/)

  end subroutine construct_flux_surfaces
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_trace_qmesh(iz, base)
  use moose_grids
  integer,     intent(in) :: iz
  type(qmesh), intent(in) :: base

  type(mmesh) :: core
  integer :: iblock, ir1


  iblock = iz / layers
  core = construct_flux_tubes(base, T(iblock)%it_base, T(iblock)%phi)

  ir1 = workspace(iz)%r_surf_pl_trans_range(1)
  workspace(iz)%r(0:ir1-1,:,:) = core%r
  workspace(iz)%z(0:ir1-1,:,:) = core%z
  workspace(iz)%b(0:ir1-1,:,:) = core%b

  end subroutine aux_trace_qmesh
  !-----------------------------------------------------------------------------
  subroutine trace_qmesh(iz, filename)
  use moose_error
  use moose_mpi
  use moose_grids
  integer,          intent(in) :: iz
  character(len=*), intent(in) :: filename

  type(qmesh) :: base
  integer :: ir1


  ir1 = workspace(iz)%r_surf_pl_trans_range(1)
  if (rank == 0) then
     print 1000, trim(filename)
     print *

     base = qmesh(filename)
     if (base%n(1) /= ir1) call ERROR("invalid radial resolution in base mesh")
     if (base%n(2) /= workspace(iz)%n(2)) call ERROR("invalid poloidal resolution in base mesh")
  endif
 1000 format(3x,"- Constructing 3-D mesh by tracing field lines from base mesh ",a)


  call base%broadcast()
  call aux_trace_qmesh(iz, base)

  end subroutine trace_qmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine mesh_generator(iz, recipe)
  !
  ! Frontend for constructing mesh in core domain of zone *iz* based on *recipe*.
  !
  use moose_error
  use moose_mpi
  use moose_quantiles
  use flare_tasks
  use flare_mmesh_parameters, only: radial_spacing
  integer,          intent(in) :: iz
  character(len=*), intent(in) :: recipe

  character(len=len(recipe)) :: cmd, dummy, err, filenames(0:blocks-1)
  class(qfunc), allocatable  :: Qr
  real(real64) :: psiN, a(0:workspace(iz)%r_surf_pl_trans_range(1)), a0
  integer :: ib, iostat, nr


  call begin_task()
  if (recipe == "") call ERROR("missing mesh generator command")
  if (rank == 0) then
     print 1000, iz, workspace(iz)%r_surf_pl_trans_range(1)
  endif
 1000 format(1x,"Constructing core domain in zone ",i0," for radial indices 0 -> ",i0/)


  read (recipe, *) cmd
  select case (cmd)
  ! construct mesh along grad psiN path
  case ("RPATH2D")
     read (recipe, *) dummy, psiN
     call construct_rpath2d(iz, psiN)


  ! construct mesh from flux surfaces
  case ("FLUX_SURFACES")
     nr = workspace(iz)%r_surf_pl_trans_range(1)
     ! explicit definition of a
     read (recipe, *, iostat=iostat) dummy, a(0:nr-1)
     ! implicit definition of a from a0 and radial_spacing(-1)
     if (iostat /= 0) then
        read (recipe, *, iostat=iostat) dummy, a0
        if (iostat /= 0) then
           call ERROR("cannot read parameter a0 for FLUX_SURFACES")
        endif
        allocate (Qr, source=generate_quantile_function(radial_spacing(-1)))
        a = a0 + (1.d0 - a0) * Qr%qquantiles(nr)
        call Qr%free();   deallocate (Qr)
     endif
     call construct_flux_surfaces(iz, a(0:nr-1))


  ! mesh is given at base, trace field lines from here
  case ("TRACE_QMESH")
     read (recipe, *) dummy, filenames
     ib = iz / layers
     call trace_qmesh(iz, filenames(ib))


  ! invalid command
  case default
     write (err, 9000) trim(cmd)
     call ERROR(err)
  end select
 9000 format("invalid mesh generator command ",a)


  call finalize_task()

  end subroutine mesh_generator
  !-----------------------------------------------------------------------------

end module flare_mmesh_core
