add_subdirectory(base)
add_subdirectory(legacy)
add_subdirectory(unstructured)
target_sources(${PROJECT_NAME} PRIVATE
    core.f90
    generator.f90
    inner_boundary.f90
    layout.f90
    mmesh.f90
    parameters.f90
    plates.f90
    read_input.f90
    utils.f90
    vacuum.f90
)
