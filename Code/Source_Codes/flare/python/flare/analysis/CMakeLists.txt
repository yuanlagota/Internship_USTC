set(fortran_analysis_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/../../../fortran/analysis)

autodoc_function(${fortran_analysis_source_dir}/melnikov_function.f90  "melnikov_function")
autodoc_function(${fortran_analysis_source_dir}/melnikov_function.f90  "melnikov_integral")
autodoc_function(${fortran_analysis_source_dir}/fluxsurf2d.f90         "equi2d_rzarray")
autodoc_function(${fortran_analysis_source_dir}/fourier_transform.f90  "toroidal_mode")
autodoc_function(${fortran_analysis_source_dir}/fourier_transform.f90  "fourier_transform")
autodoc_function(${fortran_analysis_source_dir}/fieldline.f90          "fieldline_trace")
autodoc_subroutine(${fortran_analysis_source_dir}/fluxsurf2d.f90       "fluxsurf2d_parameters")

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/${project_name}/analysis
    FILES_MATCHING PATTERN *.txt
    PATTERN "CMakeFiles" EXCLUDE
)
