add_subdirectory(_cmlib)
target_sources(${PROJECT_NAME} PRIVATE
    integral.f90
)
