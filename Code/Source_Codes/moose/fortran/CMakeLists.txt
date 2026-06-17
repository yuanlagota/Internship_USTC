# configure version module
execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION_DESCRIPTION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REGEX REPLACE "^v" "" GIT_VERSION_DESCRIPTION "${GIT_VERSION_DESCRIPTION}")
configure_file(${PROJECT_SOURCE_DIR}/config/version.f90.in version.f90 @ONLY)


# build the library
add_library(${PROJECT_NAME} SHARED)

set(fortran_legacy_compile_option
    $<$<Fortran_COMPILER_ID:GNU>:-std=legacy>
    $<$<Fortran_COMPILER_ID:Intel>:-nostand>
)

add_subdirectory(analysis)
add_subdirectory(core)
add_subdirectory(data)
add_subdirectory(geometry)
add_subdirectory(grids)
add_subdirectory(workspaces)
target_sources(${PROJECT_NAME} PRIVATE
    analysis.f90
    data.f90
    geometry.f90
    grids.f90
    workspaces.f90
    version.f90
)

target_include_directories(${PROJECT_NAME} PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${project_name}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/modules>
    ${NetCDF_Fortran_INCLUDE_DIRS}
    ${NetCDF_Fortran_INCLUDEDIR}
)

target_compile_options(${PROJECT_NAME} PRIVATE "-O2")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Og -fcheck=all -fbacktrace -Wall -ffpe-trap=invalid,zero,overflow")
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_PREPROCESS ON)

set_property(TARGET ${PROJECT_NAME} PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(${PROJECT_NAME} PUBLIC LAPACK::LAPACK PkgConfig::NetCDF_Fortran MPI::MPI_Fortran)



# install the library
install(TARGETS ${PROJECT_NAME}
    EXPORT ${project_name}-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${project_name}
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/modules/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${project_name}
)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/core/txtio.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${project_name}
)

install(EXPORT ${project_name}-targets
    FILE ${project_name}-targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
