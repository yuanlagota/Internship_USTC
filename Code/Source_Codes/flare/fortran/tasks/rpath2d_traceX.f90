subroutine flare_task_rpath2d_traceX(ix, xdir, param, t1, bounded, output)
  !
  ! Trace radial path from ix-th X-point to t1 (i.e. along :math:`\nabla \psi`).
  !
  ! **Parameters:**
  !
  ! :ix:                 X-point number.
  !
  ! :xdir:               Orientation from X-point (>>>).
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
  use flare_model, only: assert_equi2d, boundary2d, equi2d
  use flare_tasks
  use flare_rpath2d
  integer,          intent(in) :: ix, xdir
  real(real64),     intent(in) :: t1
  character(len=*), intent(in) :: param, output
  logical,          intent(in) :: bounded

  type(rpath2d_trace) :: trace
  real(real64) :: x0(2)
  integer :: iparam


  if (rank > 0) return
  call begin_task()
  call assert_equi2d("rpath2d_traceX")
  iparam = RPATH2D_PARAM(param)
  x0 = equi2d%xpoint(ix)
  print 1000, param
  print 1001, x0, xdir
 1000 format(1x,"Tracing path from X-point along grad psi direction (",a," parametrization)",/)
 1001 format(3x,"- Initial location: (",f0.3,", ",f0.3,") m, orientation: ",i0,/)


  ! construct trace
  if (bounded) then
     trace = rpath2d_traceX(ix, xdir, iparam, t1, boundary=boundary2d)
  else
     trace = rpath2d_traceX(ix, xdir, iparam, t1)
  endif


  ! save output
  call trace%savetxt(output)
  call finalize_task()

end subroutine flare_task_rpath2d_traceX
