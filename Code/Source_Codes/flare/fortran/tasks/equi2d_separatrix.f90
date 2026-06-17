subroutine flare_task_equi2d_separatrix(xpoint, output)
  !
  ! Construct separatrix for selected X-point.
  !
  ! **Parameters:**
  !
  ! :xpoint:     Separatrix is constructed for this X-point.
  !
  ! :output:     Filename output of separatrix contour.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_contours, only: xcontour, XCONTOUR_BRANCH_DIRECTION
  use flare_control
  use flare_fluxsurf2d
  use flare_tasks
  implicit none
  integer,          intent(in) :: xpoint
  character(len=*), intent(in) :: output

  type(xcontour) :: separatrix
  real(real64) :: x(2)
  integer :: i, k


  if (rank > 0) return
  ! greeting
  call begin_task()
  if (report) then
     print 1000, xpoint, trim(output)
     print *
  endif
 1000 format(1x,"Constructing separatrix for X-point ",i0," in ",a)


  ! construct separatrix and save output
  separatrix = separatrix2d(xpoint)
  call separatrix%savetxt(output)
  if (report) then
     print 1001, separatrix%x
     print *
     print 1002
     do k=0,3
        if (separatrix%iconnect(k) < 0) then
           i = 0;   if (XCONTOUR_BRANCH_DIRECTION(k) == -1) i = -1
           x = separatrix%branch(k)%point(i)
           print 1003, x, -separatrix%iconnect(k), separatrix%uconnect(k)
        endif
     enddo
  endif
 1001 format(3x,"- X-point location: (",f7.3,", ",f7.3,") m")
 1002 format(3x,"- Strike points:")
 1003 format(8x,"(",f7.3,", ",f7.3,") m on boundary ",i0," at u = ",f0.3)

  call finalize_task()

end subroutine flare_task_equi2d_separatrix
