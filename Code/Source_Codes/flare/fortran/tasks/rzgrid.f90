subroutine flare_task_rzgrid(inner_boundary1, inner_boundary2, phi, m, dr, output)
  !
  ! Construct grid in R-Z plane from poloidally closed boundaries.
  !
  ! **Parameters:**
  !
  ! :inner_boundary1:    Inner boundary contour (type(curve)).
  !
  ! :inner_boundary2:    Second inner boundary contour (optional).
  !
  ! :phi:                Toroidal angle [deg].
  !
  ! :m:                  Initial number of points in poloidal direction.
  !
  ! :dr:                 Radial resolution [m].
  !
  ! :output:             Filename for output of grid.
  !
  use iso_fortran_env
  use moose_math,       only: diff, pi
  use moose_analysis,   only: interp
  use moose_geometry,   only: polygon2d, curve, loadtxt_curve
  use moose_workspaces, only: uqwork, uqwork_contour, uqwork_layer
  use flare_boundary,   only: firstwall_rzslice
  use flare_tasks
  implicit none
  character(len=*), intent(in) :: inner_boundary1, inner_boundary2, output
  real(real64),     intent(in) :: phi, dr
  integer,          intent(in) :: m

  class(curve), allocatable :: C1, C2
  type(polygon2d) :: P
  type(uqwork) :: mesh
  real(real64) :: dp
  integer :: i, j, wn


  call begin_task()
  if (rank == 0) then
     print *, "Constructing grid in R-Z plane from poloidally closed boundaries ..."
     print *


     ! set inner boundary
     if (inner_boundary2 == "") then
        print *, "loading inner boundary ..."
        C1 = loadtxt_curve(inner_boundary1)
        mesh = uqwork_contour(C1, m, .true., iwork_nodes=1)

     else
        print *, "loading inner boundaries ..."
        C1 = loadtxt_curve(inner_boundary1)
        C2 = loadtxt_curve(inner_boundary2)
        mesh = uqwork_layer(C1, C2, m, iwork_nodes=1)
     endif
     print *


     ! set outer boundary contour
     print *, "constructing poloidally closed outer boundary from model boundary ..."
     P = firstwall_rzslice(phi / 180.d0 * pi)
     print *


     ! main loop
     print *, "adding layers along normal directoin of secants ..."
     dp = sum(norm2(diff(mesh%x(:,mesh%mnodes-m:mesh%mnodes-1), dim=2), dim=1)) / m
     i = 0
     do
        i = i + 1
        print *, i
        call mesh%add_snlayer(dr, dp)

        ! stop when all nodes are outside outer boundary
        wn = 0
        do j=mesh%lnodes,mesh%mnodes-1
           wn = wn + abs(P%winding_number(mesh%x(:,j)))
        enddo
        if (wn == 0) exit
     enddo


     print *, "saving grid in ", trim(output)
     call mesh%savetxt(output)
  endif
  call finalize_task()

end subroutine flare_task_rzgrid
