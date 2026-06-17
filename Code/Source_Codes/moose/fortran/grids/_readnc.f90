module moose_readnc_grid
  use moose_grid
  use moose_structured_grid
  implicit none

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_grid(nc) result(this)
  use moose_error
  use moose_netcdf
  use moose_mesh1d
  use moose_ugrid
  use moose_tmesh
  use moose_rmesh
  use moose_tpzmesh
  use moose_qmesh
  use moose_cmesh
  use moose_cgrid
  use moose_rmesh3d
  use moose_tpzmesh3d
  use moose_uqmesh
  type(netcdf_dataset), intent(in) :: nc
  class(grid), allocatable         :: this

  character(len=128) :: grid_type


  call nc%get_att("type", grid_type)
  select case(grid_type)
  case(TYPE_MESH1D)
     allocate (this, source = readnc_mesh1d(nc))

  case(TYPE_UGRID2D, TYPE_UGRID3D)
     allocate (this, source = readnc_ugrid(nc))

  case(TYPE_TMESH2D, TYPE_TMESH3D)
     allocate (this, source = readnc_tmesh(nc))

  case(TYPE_RMESH, TYPE_TPZMESH, TYPE_QMESH, TYPE_CMESH, TYPE_RMESH3D, TYPE_TPZMESH3D)
     allocate (this, source = readnc_structured_grid(nc))

  case(TYPE_CGRID2D, TYPE_CGRID3D)
     allocate (this, source = readnc_cgrid(nc))

  case(TYPE_UQMESH)
     allocate (this, source = readnc_uqmesh(nc))

  case default
     call ERROR("unkown grid type '"//trim(grid_type)//"'")
  end select

  end function readnc_grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function readnc_structured_grid(nc) result(this)
  use moose_error
  use moose_netcdf
  use moose_rmesh
  use moose_tpzmesh
  use moose_qmesh
  use moose_cmesh
  use moose_rmesh3d
  use moose_tpzmesh3d
  type(netcdf_dataset), intent(in) :: nc
  class(structured_grid), allocatable         :: this

  character(len=128) :: grid_type


  call nc%get_att("type", grid_type)
  select case(grid_type)
  case(TYPE_RMESH)
     allocate (this, source = readnc_rmesh(nc))

  case(TYPE_TPZMESH)
     allocate (this, source = readnc_tpzmesh(nc))

  case(TYPE_QMESH)
     allocate (this, source = readnc_qmesh(nc))

  case(TYPE_CMESH)
     allocate (this, source = readnc_cmesh(nc))

  case(TYPE_RMESH3D)
     allocate (this, source = readnc_rmesh3d(nc))

  case(TYPE_TPZMESH3D)
     allocate (this, source = readnc_tpzmesh3d(nc))

  case default
     call ERROR("unkown grid type '"//trim(grid_type)//"'")
  end select

  end function readnc_structured_grid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadnc_grid(filename) result(this)
  use moose_netcdf
  character(len=*), intent(in) :: filename
  class(grid), allocatable     :: this

  type(netcdf_dataset) :: nc


  nc = netcdf_open(filename)
  allocate (this, source = readnc_grid(nc))
  call nc%close()

  end function loadnc_grid
  !-----------------------------------------------------------------------------

end module moose_readnc_grid
