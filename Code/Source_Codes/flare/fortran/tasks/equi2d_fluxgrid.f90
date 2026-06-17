subroutine flare_task_equi2d_fluxgrid(psiN_start, psiN_end, nsteps, ntheta, grid, output)
  !
  ! Construct grid from equidistant flux coordinates.
  !
  ! **Parameters:**
  !
  ! :psiN_start:  Lower boundary [normalized poloidal flux] for radial domain of grid.
  !
  ! :psiN_end:    Upper boundary [normalized poloidal flux] for radial domain of grid.
  !
  ! :nsteps:      Number of equidistant steps between ``psiN_start`` and ``psiN_end``.
  !
  ! :ntheta:      Number of equidistant steps in poloidal direction.
  !
  ! :grid:        Filename for output of grid with R and Z coordinates.
  !
  ! :output:      Filename for output of psiN and theta values on grid nodes.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math,    only: pi2, linspace
  use moose_dataset
  use moose_grids,   only: qmesh
  use flare_control
  use flare_fluxsurf2d
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: psiN_start, psiN_end
  integer,          intent(in) :: nsteps, ntheta
  character(len=*), intent(in) :: grid, output

  real(real64), pointer :: element(:)
  real(real64)  :: psiN(0:nsteps), theta(0:ntheta), rz(2, 0:nsteps, 0:ntheta)
  type(qmesh)   :: geometry
  type(dataset) :: D
  integer :: i, k(2)


  call begin_task()
  if (rank == 0) then
     print *, "Constructing flux coordinate system for 2D equilibrium"
     print *
  endif

  ! construct flux coordinates
  psiN  = linspace(psiN_start, psiN_end, nsteps+1)
  theta = linspace(0.d0, pi2, ntheta+1)
  rz    = equi2d_rzarray(psiN, theta)

  ! write grid
  geometry = qmesh(rz(1,:,:), rz(2,:,:))
  call geometry%savetxt(grid)

  ! write Psin, theta on grid nodes
  D = dataset(2, geometry%nnodes())
  call D%set_metadata(1, "psiN",  "Normalized Poloidal Flux")
  call D%set_metadata(2, "theta", "Poloidal angle (straight field line)", "{rad}")
  do i=0,geometry%nnodes()-1
     k = geometry%node_index(i)
     element => D%element(i);   element = [psiN(k(1)), theta(k(2))]
  enddo
  call D%set_geometry(grid, output)
  call D%savetxt(output)

  call finalize_task()

end subroutine flare_task_equi2d_fluxgrid
