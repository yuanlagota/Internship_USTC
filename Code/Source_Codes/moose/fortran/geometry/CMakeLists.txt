add_subdirectory(curves)
add_subdirectory(surfaces)
target_sources(${PROJECT_NAME} PRIVATE
    polygon.f90
    polygon2d.f90
    quad.f90
    curve.f90
    surface.f90
    hypersurface.f90
    hypermesh3d.f90
    contours.f90
)
