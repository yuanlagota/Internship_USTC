subroutine flare_task_fieldline_trace(x0, x0_coordinates, x0_angular_units, direction, step_size, &
  nsteps, stop_at_boundary, trace_coordinates, angular_units, output)
  !
  ! Trace magnetic field line from x0.
  !
  ! **Parameters:**
  !
  ! :x0:                 Initial position on field line (meaning depends on `x0_coordinates`).
  !
  ! :x0_coordinates:     Coordinates for initial position `x0`:
  !
  !                      :cylindrical:  r[m], z[m], phi
  !
  !                      :cartesian:    x[m], y[m], z[m]
  !
  !                      :toroidal:     psiN, theta, phi (geometric poloidal angle)
  !
  !                      :magnetic:     psiN, theta, phi (straight field line)
  !
  ! :x0_angular_units:   Units for angular components of `x0`, if applicable (options: deg or rad).
  !
  ! :direction:          Direction of field line tracing:
  !
  !                      :fwd, forward:            in direction of toroidal field.
  !
  !                      :bwd, backward:           against direction of toroidal field.
  !
  !                      :cw, clockwise:           in clockwise direction.
  !
  !                      :ccw, counter-clockwise:  in counter-clockwise direction.
  !
  ! :step_size:          :>0: Toroidal increment [deg].
  !
  !                      :=0: Automatic set size.
  !
  ! :nsteps:             Max. number of trace steps (switched off if <= 0).
  !
  ! :stop_at_boundary:   Tracing continues beyond boundary if set to ``.false.``.
  !
  ! :trace_coordinates:  Coordinates for field line trace output (see `x0_coordinates` for options).
  !
  ! :angular_units:      Units for angular components of trace output, if applicable (options: deg or rad).
  !
  ! :output:             Filename for output of field line trace.
  !
  use iso_fortran_env
  use moose_error
  use moose_math
  use moose_units
  use flare_control
  use flare_model
  use flare_fluxsurf2d, only: fluxsurf2d
  use flare_fieldline
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: x0(3), step_size
  character(len=*), intent(in) :: x0_coordinates, x0_angular_units, direction, trace_coordinates, angular_units, output
  integer,          intent(in) :: nsteps
  logical,          intent(in) :: stop_at_boundary

  type(fluxsurf2d) :: F
  type(fieldline)  :: T
  real(real64) :: r0(3), theta, ds
  integer :: idir, n


  if (rank > 0) return
  ! 1. greeting & initial position
  call begin_task()
  if (report) then
     print *, "Tracing magnetic field line:"
     print *
     call assert_angular_units(x0_angular_units, "fieldline_trace")

     select case(x0_coordinates)
     ! cylindrical coordinates
     case(CYLINDRICAL_COORDINATES)
        r0 = x0
        select case(x0_angular_units)
        case(DEGREE)
           r0(1:2) = x0(1:2);   r0(3) = x0(3) / 180.d0 * pi
           print 1001, "r[m], z[m], phi[deg]", x0

        case(RADIAN)
           print 1001, "r[m], z[m], phi[rad]", x0
        end select


     ! Cartesian coordinates
     case(CARTESIAN_COORDINATES)
        r0 = cart_to_cyl(x0)
        print 1001, "x[m], y[m], z[m]", x0
        print 1002, r0


     ! toroidal coordinates (with respect to magnetic axis)
     case(TOROIDAL_COORDINATES)
        call assert_equi2d("fieldline_trace")
        select case(x0_angular_units)
        case(DEGREE)
           theta = x0(2) / 180.d0 * pi
           r0(3) = x0(3) / 180.d0 * pi
           print 1001, "psiN, theta[deg], phi[deg]", x0

        case(RADIAN)
           theta = x0(2)
           r0(3) = x0(3)
           print 1001, "psiN, theta[rad], phi[rad]", x0

        end select
        r0(1:2) = equi2d%rzcoords(x0(1), theta)
        print 1002, r0


     ! magnetic coordinates (straight field line)
     case("magnetic")
        call assert_equi2d("fieldline_trace")
        select case(x0_angular_units)
        case(DEGREE)
           theta = x0(2) / 180.d0 * pi
           r0(3) = x0(3) / 180.d0 * pi
           print 1001, "psiN, theta[deg], phi[deg]", x0

        case(RADIAN)
           theta = x0(2)
           r0(3) = x0(3)
           print 1001, "psiN, theta[rad], phi[rad]", x0

        end select
        F       = fluxsurf2d(equi2d%rzcoords(x0(1), 0.d0))
        r0(1:2) = F%eval(theta)
        print 1002, r0


     ! invalid input
     case default
        print 9001, x0_coordinates
        stop
     end select
  endif
  print *
 1001 format(3x,"- Initial position (",a,")"/8x,3f12.6)
 1002 format(8x,"-> in cylindrical coordinates (r[m], z[m], phi[rad]):"/8x,3f12.6)
 9001 format("ERROR in fieldline_trace: invalid x0_coordinates = '",a,"'")


  ! 2. trace parameters (steps, step size, coordinates)
  idir = make_idir(direction, bfield%equi%Bt_sign)
  if (step_size > 0.d0) then
     print 2001, abs(step_size)
  else
     print 2002
  endif
 2001 format(8x,"step size: ",f0.3," deg",/)
 2002 format(8x,"automatic step size",/)
  n = huge(1);   if (nsteps > 0) n = nsteps


  ! 3. trace field line from r0
  ds = step_size / 180.d0 * pi
  T  = fieldline(r0, idir, ds, n, stop_at_boundary, trace_coordinates, angular_units)
  if (T%at_boundary) then
     print 3001, T%trace%nelements()-1
  else
     print 3002, T%trace%nelements()-1
  endif
  call T%savetxt(output)
  call finalize_task()
 3001 format(3x,"- Field line tracing terminated at boundary after ",i0," steps")
 3002 format(3x,"- Field line tracing terminated after ",i0," steps")

end subroutine flare_task_fieldline_trace
