add_subdirectory(odeivp)
target_sources(${PROJECT_NAME} PRIVATE
    fieldline.f90
    qmesh_generator.f90
    uqwork.f90
    bfit.f90
)
