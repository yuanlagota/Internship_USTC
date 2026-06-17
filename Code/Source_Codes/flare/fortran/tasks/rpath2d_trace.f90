subroutine flare_task_rpath2d_trace(x0, param, t1, bounded, output)
  !
  ! Trace radial path from x0 (i.e. along :math:`\nabla \psi`).
  !
  ! **Parameters:**
  !
  ! :x0:                 Initial position (R[m], Z[m]).
  !
  ! :param:              Parametrization of trace:
  !
  !                      :arclength:    arc length [m].
  !
  !                      :psiN:         normalized poloidal flux.
  !
  ! :t1:                 Either max. arc length or destination psiN.
  !
  ! :bounded:            Truncate trace at model boundary if set.
  !
  ! :output:             Filename for output for trace.
  !
  use iso_fortran_env
  use flare_model, only: assert_equi2d, boundary2d
  use flare_tasks
  use flare_rpath2d
  real(real64),     intent(in) :: x0(2), t1
  character(len=*), intent(in) :: param, output
  logical,          intent(in) :: bounded

  type(rpath2d_trace) :: trace
  integer :: iparam


  if (rank > 0) return
  call begin_task()
  call assert_equi2d("rpath2d_trace")
  iparam = RPATH2D_PARAM(param)
  print 1000, param
  print 1001, x0
 1000 format(1x,"Tracing path along grad psi direction (",a," parametrization)",/)
 1001 format(3x,"- Initial location: (",f0.3,", ",f0.3,") m",/)


  ! construct trace
  if (bounded) then
     trace = rpath2d_trace(x0, iparam, t1, boundary=boundary2d)
  else
     trace = rpath2d_trace(x0, iparam, t1)
  endif


  ! save output
  call trace%savetxt(output)
  call finalize_task()

end subroutine flare_task_rpath2d_trace
