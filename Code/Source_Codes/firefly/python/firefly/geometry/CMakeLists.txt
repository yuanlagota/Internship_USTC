set(fortran_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/../../../fortran)

autodoc_subroutine(${fortran_source_dir}/geometry.f90  "init_workspace")
autodoc_subroutine(${fortran_source_dir}/geometry.f90  "set_pfc")

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/${PROJECT_NAME}/geometry
    FILES_MATCHING PATTERN *.txt
    PATTERN "CMakeFiles" EXCLUDE
)
