add_subdirectory(math)
add_subdirectory(utils)
target_sources(${PROJECT_NAME} PRIVATE
    aux.c
    kinds.f90
    error.f90
    mpi.f90
    netcdf.f90
    algorithms.f90
    units.f90
    utils.f90
    math.f90
    dict.f90
    configparser.f90
    txtio.f90
    table.f90
    ilist.f90
    rlist.f90
)

# add include directory for txtio.h
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
