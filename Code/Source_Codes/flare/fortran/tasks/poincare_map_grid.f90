subroutine flare_task_poincare_map_grid(grid, direction, phi_section, nsym, npoints, nsections, bounded, output)
  !
  ! Create Poincare maps for initial points from ``grid``.
  !
  ! **Parameters:**
  !
  ! :grid:         Filename for grid with initial points.
  !
  ! :direction:    Direction of field line tracing:
  !
  !                :fwd, forward:            in direction of toroidal field.
  !
  !                :bwd, backward:           against direction of toroidal field.
  !
  !                :cw, clockwise:           in clockwise direction.
  !
  !                :ccw, counter-clockwise:  in counter-clockwise direction.
  !
  ! :phi_section:  Toroidal position of Poincare section [deg].
  !
  ! :nsym:         Toroidal symmetry of Poincare section.
  !
  ! :npoints:      Max. number of return points in Poincare section for each starting point.
  !
  ! :nsections:    Number of Poincare sections with toroidal positions at :math:`\varphi_0 \, \pm \, i \, \cdot \, \frac{360 \, \deg}{n_{sym} \, n_{sections}}, \qquad i = 0 \ldots n_{sections} - 1`.
  !
  ! :bounded:      Field line tracing terminates at boundary.
  !
  ! :output:       Filename for output of Poincare maps.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_grids, only: r3grid
  implicit none
  character(len=*), intent(in) :: grid, direction, output
  real(real64),     intent(in) :: phi_section
  integer,          intent(in) :: nsym, npoints, nsections
  logical,          intent(in) :: bounded

  type(r3grid) :: G
  character(len=256) :: label


  if (rank == 0) then
     G = r3grid(grid)
     write (label, 1003) G%nnodes(), grid
  endif
  call G%broadcast()
 1003 format(i0," initial points from ",a)


  call flare_task_poincare_map_grid_implementation(G, label, direction, phi_section, nsym, npoints, nsections, bounded, output)

end subroutine flare_task_poincare_map_grid





subroutine flare_task_poincare_map_grid_implementation(G, label, direction, phi0, nsym, npoints, nsections, bounded, output)
  use iso_fortran_env
  use moose_mpi
  use moose_utils, only: str
  use moose_math,  only: pi, linspace
  use moose_grids, only: r3grid
  use flare_control
  use flare_model, only: bfield, assert_equi2d, equi2d
  use flare_fieldline
  use flare_poincare_map
  use flare_tasks
  implicit none
  type(r3grid),     intent(in) :: G
  character(len=*), intent(in) :: label, direction, output
  real(real64),     intent(in) :: phi0
  integer,          intent(in) :: nsym, npoints, nsections
  logical,          intent(in) :: bounded

  type(poincare_map), allocatable :: P(:,:)
  character(len=120)    :: filename, suffix
  type(fdriver), target :: F
  real(real64)     :: x0(3), phi0_rad
  character(len=1) :: c
  integer :: i, idir, imod, iu, j, k, n, m


  ! greeting & task parameters
  call begin_task()
  if (rank == 0) then
     print *, "Generating Poincare maps ..."
     print *

     idir = bfield%equi%Bt_sign * make_idir(direction, bfield%equi%Bt_sign)
     print 1001, phi0
     print 1002, nsym
     print *
     print 1003, trim(label)
     print *
  endif
  call proc(0)%broadcast(idir)
  phi0_rad = phi0 / 180.d0 * pi
 1001 format(8x,"reference position: ",f0.3," deg")
 1002 format(8x,"toroidal symmetry:  ",i0)
 1003 format(3x,"- ",a)


  ! compute Poincare maps
  F = fdriver(stop_at_boundary=bounded)
  m = G%nnodes()
  allocate (P(0:nsections-1,0:m-1))
  do i=rank,m-1,nproc
     x0   = G%node(i)
     P(:,i) = poincare_maps(x0, idir, phi0_rad, nsym, npoints, nsections, fdriver=F)

     ! display summary
     n = P(0,i)%points%nelements()
     c = "=";   if (n == npoints) c = ">"
     if (P(0,i)%istat > 0) c = str(P(0,i)%istat)
     print 3000, i, c, F%arclength, n
  enddo
 3000 format(3x,i8,":",3x,"Lc ",a,3x,f9.2," m",4x,"with ",i4," points")


  ! collect results & save output
  do j=0,nsections-1
     filename = output
     if (nsections > 1) then
        suffix = ""
        k = index(output, ".", back=.true.)
        if (k > 0) then
           filename = output(1:k-1)
           suffix   = output(k:)
        endif
        write (filename, 4000) trim(filename), j, trim(suffix)
     endif
     if (rank == 0) then
        open  (newunit=iu, file=filename)
        write (iu, 4001)
     endif

     do i=0,m-1
        imod = mod(i,nproc)
        ! send results to first rank
        if (imod /= 0) then
           if (imod == rank) call P(j,i)%send(0)
           if (rank == 0) P(j,i) = recv_poincare_map(imod)
        endif

        ! write output
        if (rank == 0) write (iu, '(dt)') P(j,i)
     enddo
     if (rank == 0) close(iu)
  enddo
  call finalize_task()
  deallocate (P)
 4000 format(a,"_",i0,a)
 4001 format("# TYPE poincare_map")

end subroutine flare_task_poincare_map_grid_implementation
