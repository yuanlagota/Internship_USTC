!===============================================================================
! ODE stepper for embedded Runga-Kutta methods
!
! References:
!             A family of embedded Runge-Kutta formulae, Dormand & Prince
!             Journal of Computational and Applied Mathematics, Vol. 6, 19-26, 1980.
!
!             High order embedded Runge-Kutta formulae, Prince & Dormand,
!             Journal of Computational and Applied Mathematics, Vol. 7, 67-75, 1981
!
!             E. Hairer, S.P. Norsett, G. Wanner, "Solving ordinary differential equations I,
!             Nonstiff Problems", 2nd revised edition, Springer, 2000
!===============================================================================
module moose_odeivp_rk_plugin
  use iso_fortran_env
  use moose_error
  use moose_odeivp_system
  use moose_odeivp_stepper
  implicit none
  private


  integer, parameter, public :: &
     STEP_INCREASE =  1, &
     STEP_DECREASE = -1


  ! Runge-Kutta-Fehlberg method (rkf45)
  real(real64), parameter :: &
     rkf45_c(2:6)  = [1.d0/4.d0, 3.d0/8.d0, 12.d0/13.d0, 1.d0, 1.d0/2.d0], &
     rkf45_a2(1)   = [1.d0/4.d0], &
     rkf45_a3(2)   = [3.d0/32.d0, 9.d0/32.d0], &
     rkf45_a4(3)   = [1932.d0/2197.d0, -7200.d0/2197.d0, 7296.d0/2197.d0], &
     rkf45_a5(4)   = [439.d0/216.d0, -8.d0, 3680.d0/513.d0, -845.d0/4104.d0], &
     rkf45_a6(5)   = [-8/27.d0, 2.d0, -3544.d0/2565.d0, 1859.d0/4104.d0, -11.d0/40.d0], &
     rkf45_bhat(6) = [16.d0/135.d0, 0.d0, 6656.d0/12825.d0, 28561.d0/56430.d0, -9.d0/50.d0, 2.d0/55.d0], &
     rkf45_b(6)    = [25.d0/216.d0, 0.d0, 1408.d0/2565.d0, 2197.d0/4104.d0, -1.d0/5.d0, 0.d0]


  ! RK5(4)7M (dopr5)
  real(real64), parameter :: &
     rk547m_c(2:7)  = [1.d0/5, 3.d0/10, 4.d0/5, 8.d0/9, 1.d0, 1.d0], &
     rk547m_a2(1)   = [1.d0/5], &
     rk547m_a3(2)   = [3.d0/40, 9.d0/40], &
     rk547m_a4(3)   = [44.d0/45, -56.d0/15, 32.d0/9], &
     rk547m_a5(4)   = [19372.d0/6561, -25360.d0/2187, 64448.d0/6561, -212.d0/729], &
     rk547m_a6(5)   = [9017.d0/3168, -355.d0/33, 46732.d0/5247, 49.d0/176, -5103.d0/18656], &
     rk547m_a7(6)   = [35.d0/384, 0.d0, 500.d0/1113, 125.d0/192, -2187.d0/6784, 11.d0/84], &
     rk547m_bhat(7) = [rk547m_a7, 0.d0], &
     rk547m_b(7)    = [5179.d0/57600, 0.d0, 7571.d0/16695, 393.d0/640, -92097.d0/339200, 187.d0/2100, 1.d0/40]


  ! RK5(4)7C
  real(real64), parameter :: &
     rk547c_c(2:7)  = [1.d0/5, 3.d0/10, 6.d0/13, 2.d0/3, 1.d0, 1.d0], &
     rk547c_a2(1)   = [1.d0/5], &
     rk547c_a3(2)   = [3.d0/40, 9.d0/40], &
     rk547c_a4(3)   = [264.d0/2197, -90.d0/2197, 840.d0/2197], &
     rk547c_a5(4)   = [932.d0/3645, -14.d0/27, 3256.d0/5103, 7436.d0/25515], &
     rk547c_a6(5)   = [-367.d0/513, 30.d0/19, 9940.d0/5643, -29575.d0/8208, 6615.d0/3344], &
     rk547c_a7(6)   = [35.d0/432, 0.d0, 8500.d0/14553, -28561.d0/84672, 405.d0/704, 19.d0/196], &
     rk547c_bhat(7) = [rk547c_a7, 0.d0], &
     rk547c_b(7)    = [11.d0/108, 0.d0, 6250.d0/14553, -2197.d0/21168, 81.d0/176, 171.d0/1960, 1.d0/40]


  ! RK5(4)7S
  real(real64), parameter :: &
     rk547s_c(2:7)  = [2.d0/9, 1.d0/3, 5.d0/9, 2.d0/3, 1.d0, 1.d0], &
     rk547s_a2(1)   = [2.d0/9], &
     rk547s_a3(2)   = [1.d0/12, 1.d0/4], &
     rk547s_a4(3)   = [55.d0/324, -25.d0/108, 50.d0/81], &
     rk547s_a5(4)   = [83.d0/330, -13.d0/22, 61.d0/66, 9.d0/110], &
     rk547s_a6(5)   = [-19.d0/28, 9.d0/4, 1.d0/7, -27.d0/7, 22.d0/7], &
     rk547s_a7(6)   = [19.d0/200, 0.d0, 3.d0/5, -243.d0/400, 33.d0/40, 7.d0/80], &
     rk547s_bhat(7) = [rk547s_a7, 0.d0], &
     rk547s_b(7)    = [431.d0/5000, 0.d0, 333.d0/500, -7857.d0/10000, 957.d0/1000, 193.d0/2000, -1.d0/50]


  ! RK6(5)8M
  real(real64), parameter :: &
     rk658m_c(2:8)  = [0.1d0, 2.d0/9, 3.d0/7, 3.d0/5, 4.d0/5, 1.d0, 1.d0], &
     rk658m_a2(1)   = [0.1d0], &
     rk658m_a3(2)   = [-2.d0/81, 20.d0/81], &
     rk658m_a4(3)   = [615.d0/1372, -270.d0/343, 1053.d0/1372], &
     rk658m_a5(4)   = [3243.d0/5500, -54.d0/55, 50949.d0/71500, 4998.d0/17875], &
     rk658m_a6(5)   = [-26492.d0/37125, 72.d0/55, 2808.d0/23375, -24206.d0/37125, 338.d0/459], &
     rk658m_a7(6)   = [5561.d0/2376, -35.d0/11, -24117.d0/31603, 899983.d0/200772, -5225.d0/1836, 3925.d0/4056], &
     rk658m_a8(7)   = [465467.d0/266112, -2945.d0/1232, -5610201.d0/14158144, 10513573.d0/3212352, -424325.d0/205632, 376225.d0/454272, 0.d0], &
     rk658m_bhat(8) = [61.d0/864, 0.d0, 98415.d0/321776, 16807.d0/146016, 1375.d0/7344, 1375.d0/5408, -37.d0/1120, 0.1d0], &
     rk658m_b(8)    = [821.d0/10800, 0.d0, 19683.d0/71825, 175273.d0/912600, 395.d0/3672, 785.d0/2704, 3.d0/50, 0.d0]


  ! RK6(5)8C
  real(real64), parameter :: &
     rk658c_c(2:8)  = [0.1d0, 1.d0/6, 2.d0/9, 3.d0/5, 4.d0/5, 1.d0, 1.d0], &
     rk658c_a2(1)   = [0.1d0], &
     rk658c_a3(2)   = [1.d0/36, 5.d0/36], &
     rk658c_a4(3)   = [10.d0/243, 20.d0/243, 8.d0/81], &
     rk658c_a5(4)   = [4047.d0/5500, -18.d0/55, -4212.d0/1375, 17901.d0/5500], &
     rk658c_a6(5)   = [-5587.d0/4125, 24.d0/55, 9576.d0/1375, -140049.d0/23375, 38.d0/51], &
     rk658c_a7(6)   = [12961.d0/2376, -35.d0/33, -160845.d0/5434, 1067565.d0/38896, -103375.d0/47736, 32875.d0/35568], &
     rk658c_a8(7)   = [702799.d0/199584, -1865.d0/2772, -2891375.d0/152152, 19332955.d0/1089088, -5356375.d0/4009824, 2207875.d0/2987712, 0.d0], &
     rk658c_bhat(8) = [1.d0/12, 0.d0, -216.d0/1235, 6561.d0/12376, 1375.d0/5304, 1375.d0/5928, -5.d0/168, 0.1d0], &
     rk658c_b(8)    = [163.d0/1440, 0.d0, -2628.d0/6175, 13851.d0/17680, 1525.d0/7956, 6575.d0/23712, 3.d0/50, 0.d0]


  ! RK6(5)8s
  real(real64), parameter :: &
     rk658s_c(2:8)  = [1.d0/4, 3.d0/10, 6.d0/7, 3.d0/5, 4.d0/5, 1.d0, 1.d0], &
     rk658s_a2(1)   = [1.d0/4], &
     rk658s_a3(2)   = [3.d0/25, 9.d0/50], &
     rk658s_a4(3)   = [102.d0/343, -1368.d0/343, 1560.d0/343], &
     rk658s_a5(4)   = [-3.d0/100, 36.d0/25, -12.d0/13, 147.d0/1300], &
     rk658s_a6(5)   = [37.d0/225, -48.d0/25, 872.d0/351, 49.d0/1053, 2.d0/81], &
     rk658s_a7(6)   = [11.d0/648, 14.d0/3, -10193.d0/2106, -30331.d0/50544, 1025.d0/1944, 59.d0/48], &
     rk658s_a8(7)   = [796.d0/1701, -352.d0/63, 134093.d0/22113, -78281.d0/75816, -9425.d0/20412, 781.d0/504, 0.d0], &
     rk658s_bhat(8) = [29.d0/324, 0.d0, 3400.d0/7371, -16807.d0/25272, -125.d0/1944, 25.d0/24, 1.d0/84, 1.d0/8], &
     rk658s_b(8)    = [2041.d0/21600, 0.d0, 748.d0/1755, -2401.d0/46800, 11.d0/108, 59.d0/160, 3.d0/50, 0.d0]


  ! RK8(7)13M
  real(real64), parameter :: &
     rk8713m_c(2:13)  = [1.d0/18, 1.d0/12, 1.d0/8, 5.d0/16, 3.d0/8, 59.d0/400, 93.d0/200, 5490023248.d0/9719169821.d0, 13.d0/20, 1201146811.d0/1299019798.d0, 1.d0, 1.d0], &
     rk8713m_a2(1)    = [1.d0/18], &
     rk8713m_a3(2)    = [1.d0/48, 1.d0/16], &
     rk8713m_a4(3)    = [1.d0/32, 0.d0, 3.d0/32], &
     rk8713m_a5(4)    = [5.d0/16, 0.d0, -75.d0/64, 75.d0/64], &
     rk8713m_a6(5)    = [3.d0/80, 0.d0, 0.d0, 3.d0/16, 3.d0/20], &
     rk8713m_a7(6)    = [29443841.d0/614563906.d0, 0.d0, 0.d0, 77736538.d0/692538347.d0, -28693883.d0/1125000000.d0, 23124283.d0/1800000000.d0], &
     rk8713m_a8(7)    = [16016141.d0 / 946692911.d0, &
                         0.d0, &
                         0.d0, &
                         61564180.d0 / 158732637.d0, &
                         22789713.d0 / 633445777.d0, &
                         545815736.d0 / 2771057229.d0, &
                         -180193667.d0 / 1043307555.d0], &
     rk8713m_a9(8)    = [39632708.d0 / 573591083.d0, &
                         0.d0, &
                         0.d0, &
                         -433636366.d0 / 683701615.d0, &
                         -421739975.d0 / 2616292301.d0, &
                         100302831.d0 / 723423059.d0, &
                         790204164.d0 / 839813087.d0, &
                         800635310.d0 / 3783071287.d0], &
     rk8713m_a10(9)   = [246121993.d0 / 1340847787.d0, &
                         0.d0, &
                         0.d0, &
                         -37695042795.d0 / 15268766246.d0, &
                         -309121744.d0 / 1061227803.d0, &
                         -12992083.d0 / 490766935.d0, &
                         6005943493.d0 / 2108947869.d0, &
                         393006217.d0 / 1396673457.d0, &
                         123872331.d0 / 1001029789.d0], &
     rk8713m_a11(10)  = [-1028468189.d0 / 846180014.d0, &
                         0.d0, &
                         0.d0, &
                         8478235783.d0 / 508512852.d0, &
                         1311729495.d0 / 1432422823.d0, &
                         -10304129995.d0 / 1701304382.d0, &
                         -48777925059.d0 / 3047939560.d0, &
                         15336726248.d0 / 1032824649.d0, &
                         -45442868181.d0 / 3398467696.d0, &
                         3065993473.d0 / 597172653.d0], &
     rk8713m_a12(11)  = [185892177.d0 / 718116043.d0, &
                         0.d0, &
                         0.d0, &
                         -3185094517.d0 / 667107341.d0, &
                         -477755414.d0 / 1098053517.d0, &
                         -703635378.d0 / 230739211.d0, &
                         5731566787.d0 / 1027545527.d0, &
                         5232866602.d0 / 850066563.d0, &
                         -4093664535.d0 / 808688257.d0, &
                         3962137247.d0 / 1805957418.d0, &
                         65686358.d0 / 487910083.d0], &
     rk8713m_a13(12)  = [403863854.d0 / 491063109.d0, &
                         0.d0, &
                         0.d0, &
                         -5068492393.d0 / 434740067.d0, &
                         -411421997.d0 / 543043805.d0, &
                         652783627.d0 / 914296604.d0, &
                         11173962825.d0 / 925320556.d0, &
                         -13158990841.d0 / 6184727034.d0, &
                         3936647629.d0 / 1978049680.d0, &
                         -160528059.d0 / 685178525.d0, &
                         248638103.d0 / 1413531060.d0, &
                         0.d0], &
     rk8713m_bhat(13) = [14005451.d0 / 335480064.d0, &
                         0.d0, &
                         0.d0, &
                         0.d0, &
                         0.d0, &
                         -59238493.d0 / 1068277825.d0, &
                         181606767.d0 / 758867731.d0, &
                         561292985.d0 / 797845732.d0, &
                         -1041891430.d0 / 1371343529.d0, &
                         760417239.d0 / 1151165299.d0, &
                         118820643.d0 / 751138087.d0, &
                         -528747749.d0 / 2220607170.d0, &
                         1.d0 / 4.d0], &
     rk8713m_b(13)    = [13451932.d0 / 455176623.d0, &
                         0.d0, &
                         0.d0, &
                         0.d0, &
                         0.d0, &
                         -808719846.d0 / 976000145.d0, &
                         1757004468.d0 / 5645159321.d0, &
                         656045339.d0 / 265891186.d0, &
                         -3867574721.d0 / 1518517206.d0, &
                         465885868.d0 / 322736535.d0, &
                         53011238.d0 / 667516719.d0, &
                         2.d0 / 45.d0, &
                         0.d0]



  ! stepper program for embedded Runge-Kutta methods ...........................
  type, extends(odeivp_stepper), public :: rk_plugin
     ! step method
     real(real64), allocatable :: a(:,:), b(:), c(:), e(:)
     integer :: nstage
     logical :: fsal

     ! controller parameters
     real(real64) :: ay, ayprime

     contains
     procedure :: free, errlevel, control, step, try_step
  end type rk_plugin


  interface rk_plugin
     procedure :: init_rk_plugin
  end interface rk_plugin
  ! rk_plugin ..................................................................


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function init_rk_plugin(system, method, hstart, epsabs, epsrel, hmin, hmax) result(this)
  use moose_utils, only: str
  class(odeivp_system), target, intent(in) :: system
  character(len=*),             intent(in) :: method
  real(real64),                 intent(in) :: hstart, epsabs, epsrel
  real(real64),                 intent(in), optional :: hmin, hmax
  type(rk_plugin)                          :: this


  ! local extrapolation -> use higher order method for numerical result, even though
  ! the error estimate is for the lower order method
  select case(method)
  case('rkf45')
     this%norder = 5
     this%nstage = 6
     this%fsal = .false.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rkf45_bhat)
     allocate (this%e, source=rkf45_bhat-rkf45_b)
     allocate (this%c(2:6), source=rkf45_c)
     this%a(2,1:1) = rkf45_a2
     this%a(3,1:2) = rkf45_a3
     this%a(4,1:3) = rkf45_a4
     this%a(5,1:4) = rkf45_a5
     this%a(6,1:5) = rkf45_a6


  case('rk5(4)7m', 'dopr5')
     this%norder = 5
     this%nstage = 7
     this%fsal = .true.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rk547m_bhat)
     allocate (this%e, source=rk547m_bhat-rk547m_b)
     allocate (this%c(2:7), source=rk547m_c)
     this%a(2,1:1) = rk547m_a2
     this%a(3,1:2) = rk547m_a3
     this%a(4,1:3) = rk547m_a4
     this%a(5,1:4) = rk547m_a5
     this%a(6,1:5) = rk547m_a6
     this%a(7,1:6) = rk547m_a7


  case('rk5(4)7c')
     this%norder = 5
     this%nstage = 7
     this%fsal = .true.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rk547c_bhat)
     allocate (this%e, source=rk547c_bhat-rk547c_b)
     allocate (this%c(2:7), source=rk547c_c)
     this%a(2,1:1) = rk547c_a2
     this%a(3,1:2) = rk547c_a3
     this%a(4,1:3) = rk547c_a4
     this%a(5,1:4) = rk547c_a5
     this%a(6,1:5) = rk547c_a6
     this%a(7,1:6) = rk547c_a7


  case('rk5(4)7s')
     this%norder = 5
     this%nstage = 7
     this%fsal = .true.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rk547s_bhat)
     allocate (this%e, source=rk547s_bhat-rk547s_b)
     allocate (this%c(2:7), source=rk547s_c)
     this%a(2,1:1) = rk547s_a2
     this%a(3,1:2) = rk547s_a3
     this%a(4,1:3) = rk547s_a4
     this%a(5,1:4) = rk547s_a5
     this%a(6,1:5) = rk547s_a6
     this%a(7,1:6) = rk547s_a7


  case('rk6(5)8m', 'dopr6')
     this%norder = 6
     this%nstage = 8
     this%fsal = .false.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rk658m_bhat)
     allocate (this%e, source=rk658m_bhat-rk658m_b)
     allocate (this%c(2:8), source=rk658m_c)
     this%a(2,1:1) = rk658m_a2
     this%a(3,1:2) = rk658m_a3
     this%a(4,1:3) = rk658m_a4
     this%a(5,1:4) = rk658m_a5
     this%a(6,1:5) = rk658m_a6
     this%a(7,1:6) = rk658m_a7
     this%a(8,1:7) = rk658m_a8


  case('rk6(5)8c')
     this%norder = 6
     this%nstage = 8
     this%fsal = .false.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rk658c_bhat)
     allocate (this%e, source=rk658c_bhat-rk658c_b)
     allocate (this%c(2:8), source=rk658c_c)
     this%a(2,1:1) = rk658c_a2
     this%a(3,1:2) = rk658c_a3
     this%a(4,1:3) = rk658c_a4
     this%a(5,1:4) = rk658c_a5
     this%a(6,1:5) = rk658c_a6
     this%a(7,1:6) = rk658c_a7
     this%a(8,1:7) = rk658c_a8


  case('rk6(5)8s')
     this%norder = 6
     this%nstage = 8
     this%fsal = .false.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rk658s_bhat)
     allocate (this%e, source=rk658s_bhat-rk658s_b)
     allocate (this%c(2:8), source=rk658s_c)
     this%a(2,1:1) = rk658s_a2
     this%a(3,1:2) = rk658s_a3
     this%a(4,1:3) = rk658s_a4
     this%a(5,1:4) = rk658s_a5
     this%a(6,1:5) = rk658s_a6
     this%a(7,1:6) = rk658s_a7
     this%a(8,1:7) = rk658s_a8


  case('rk8(7)13m', 'dopr8')
     this%norder = 8
     this%nstage = 13
     this%fsal = .false.

     allocate (this%a(2:this%nstage,this%nstage-1), source=0.d0)
     allocate (this%b, source=rk8713m_bhat)
     allocate (this%e, source=rk8713m_bhat-rk8713m_b)
     allocate (this%c(2:13), source=rk8713m_c)
     this%a(2,1:1)   = rk8713m_a2
     this%a(3,1:2)   = rk8713m_a3
     this%a(4,1:3)   = rk8713m_a4
     this%a(5,1:4)   = rk8713m_a5
     this%a(6,1:5)   = rk8713m_a6
     this%a(7,1:6)   = rk8713m_a7
     this%a(8,1:7)   = rk8713m_a8
     this%a(9,1:8)   = rk8713m_a9
     this%a(10,1:9)  = rk8713m_a10
     this%a(11,1:10) = rk8713m_a11
     this%a(12,1:11) = rk8713m_a12
     this%a(13,1:12) = rk8713m_a13


  case default
     call ERROR("invalid method '"//trim(method)//"'", "init_rk_plugin")

  end select


  this%ay = 1.d0
  this%ayprime = 0.d0
  call aux_init_odeivp_stepper(this, system, hstart, epsabs, epsrel, hmin, hmax)
  if (this%hmax == 0.d0) this%hmax = huge(1.d0)

  end function init_rk_plugin
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(rk_plugin), intent(inout) :: this


  deallocate (this%a, this%b, this%c, this%e)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function errlevel(this, y, yprime, h, errlev) result(istat)
  !
  ! calculate desired error level *errlev*.
  !
  class(rk_plugin), intent(in   ) :: this
  real(real64),     intent(in   ) :: y(:), yprime(size(y)), h
  real(real64),     intent(  out) :: errlev(size(y))
  integer                         :: istat


  istat  = SUCCESS
  errlev = this%epsabs + this%epsrel * (this%ay * abs(y) + this%ayprime * abs(h * yprime))
  if (any(errlev <= 0.d0)) istat = SANITY_ERROR

  end function errlevel
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function control(this, y, yerr, yprime, h) result(istat)
  !
  ! adjust the step-size h based on current values y, yerr and yprim
  !
  class(rk_plugin), intent(in   ) :: this
  real(real64),     intent(in   ) :: y(this%ndim), yerr(this%ndim), yprime(this%ndim)
  real(real64),     intent(inout) :: h
  integer                         :: istat

  real(real64), parameter :: &
     safety = 0.9d0, &  ! safety factor
     fshrnk = 0.2d0, &  ! limit for reducing step size
     fgrow  = 5.d0      ! limit for increasing step size

  real(real64) :: scal(size(y)), errmax


  ! evaluate tolerance
  istat = this%errlevel(y, yprime, h, scal);   if (istat /= SUCCESS) return
  errmax = maxval(abs(yerr) / scal)


  ! observed error exceeds desired error by more than 10% in at least one component
  if (errmax > 1.1d0) then
     ! decrease step-size (but no more than by fshrnk)
     h = h * max(safety * errmax**(-1.d0/this%norder), fshrnk)
     istat = STEP_DECREASE


  ! observed error is less than 50% of the desired error level
  else if (errmax < 0.5d0) then
     ! increase step-size (but no more than by factor fgrow)
     h = h * max(1.d0, min(safety * errmax**(-1.d0/(this%norder+1)), fgrow))
     istat = STEP_INCREASE
  endif

  end function control
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function step(this, t, t1, y, yprime) result(istat)
  use moose_math, only: sign_test
  class(rk_plugin), intent(inout) :: this
  real(real64),     intent(inout) :: t, y(this%ndim), yprime(this%ndim)
  real(real64),     intent(in   ) :: t1
  integer                         :: istat

  real(real64) :: dt, htry, hnext, hstart, t0, tnext, y0(this%ndim), yprime0(this%ndim), yerr(this%ndim)
  logical :: final_step


  ! initialize first call
  dt = t1 - t
  if (this%nsteps == 0) then
     hstart = min(this%hstart, this%hmax)
     this%h = sign(hstart, dt)

  ! check integration direction
  elseif (sign_test(dt, this%h) < 0) then
     istat = 3
     return
  endif
  t0 = t
  y0 = y
  yprime0 = yprime


  htry = this%h
  try_step: do
     final_step = .false.
     ! update step size htry if this is the final step
     if (sign_test(dt, htry-dt) >= 0) then
        htry = dt
        final_step = .true.
     endif

     ! execute step
     istat = this%try_step(t0, htry, y, yerr, yprime)
     ! return if stepper indicates a non-recoverable error
     if (istat == USER_FUNCTION_ERROR) return
     ! reduce step size after recoverable user function error
     if (istat /= SUCCESS) then
        hnext = htry / 2
        ! verify finite - but decreased - step size
        tnext = t + hnext
        if (abs(hnext) < abs(htry)  .and.  t /= tnext  .and.  abs(hnext) > this%hmin) then
           ! undo step and try again
           y = y0
           yprime = yprime0
           htry = hnext
           this%nfailed_steps = this%nfailed_steps + 1
           cycle
        else
           ! return and notify user of the failed step size
           this%h = htry
           t = t0
           return
        endif
     endif

     ! update current t
     t = t0 + htry
     if (final_step) t = t1

     ! check error and attempt to adjust the step
     hnext = htry
     istat = this%control(y, yerr, yprime, hnext)
     if (istat == STEP_DECREASE) then
        ! verify finite - but decreased - step size
        tnext = t + hnext
        if (abs(hnext) < abs(htry)  .and.  t /= tnext  .and.  abs(hnext) > this%hmin) then
           ! undo step and try again
           y = y0
           yprime = yprime0
           htry = hnext
           this%nfailed_steps = this%nfailed_steps + 1
           cycle
        else
           ! notify user of the failed step size
           this%h = hnext
           istat = NO_PROGRESS_ERROR
           return
        endif
     elseif (istat == STEP_INCREASE  .and.  abs(hnext) > this%hmax) then
        hnext = sign(this%hmax, hnext)
     endif

     ! suggest new step size for next time step
     if (.not.final_step) this%h = hnext
     istat = SUCCESS
     this%nsteps = this%nsteps + 1
     exit try_step
  enddo try_step

  end function step
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function try_step(this, t, h, y, yerr, yprime) result(istat)
  class(rk_plugin), intent(inout) :: this
  real(real64),     intent(in   ) :: t, h
  real(real64),     intent(inout) :: y(this%ndim), yprime(this%ndim)
  real(real64),     intent(  out) :: yerr(this%ndim)
  integer                         :: istat

  real(real64) :: ytmp(this%ndim), d(this%ndim), k(this%ndim, this%nstage)
  integer :: i, j


  ! k1
  k(:,1) = yprime


  ! k(i)
  do i=2,this%nstage
     d = 0.d0
     do j=1,i-1
        d = d + this%a(i,j) * k(:,j)
     enddo
     ytmp = y + h*d

     istat = this%system%eval(t + this%c(i)*h, ytmp, k(:,i))
     if (istat /= SUCCESS) return
  enddo


  ! final sum
  if (this%fsal) then
     y = ytmp
     yprime = k(:,this%nstage)

  else
     d = 0.d0
     do j=1,this%nstage
        d = d + this%b(j) * k(:,j)
     enddo
     y = y + h*d

     ! yprime output for next step
     istat = this%system%eval(t + h, y, yprime)
     if (istat /= SUCCESS) return
  endif


  ! error estimate
  yerr = 0.d0
  do i=1,this%nstage
     yerr = yerr + this%e(i) * k(:,i)
  enddo
  yerr = yerr * h

  end function try_step
  !-----------------------------------------------------------------------------

end module moose_odeivp_rk_plugin
