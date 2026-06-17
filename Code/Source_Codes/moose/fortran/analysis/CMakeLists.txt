add_subdirectory(mfuncs)
add_subdirectory(ufuncs)
target_sources(${PROJECT_NAME} PRIVATE
    linalg.f90
    ufunc.f90
    mfunc.f90
    quantiles.f90
)
