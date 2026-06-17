set(fortran_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/../../../fortran)

autodoc_subroutine(${fortran_source_dir}/tasks.f90  "connection_length")
autodoc_subroutine(${fortran_source_dir}/tasks.f90  "strike_point_density")
autodoc_subroutine(${fortran_source_dir}/tasks.f90  "heat_load_proxy")

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/${PROJECT_NAME}/tasks
    FILES_MATCHING PATTERN *.txt
    PATTERN "CMakeFiles" EXCLUDE
)
