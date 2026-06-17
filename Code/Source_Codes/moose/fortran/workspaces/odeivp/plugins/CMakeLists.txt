set(dlsode_source_file _dlsode.f)

target_sources(${PROJECT_NAME} PRIVATE
    ${dlsode_source_file}
    rk.f90
    lsode.f90
)

set_source_files_properties(${dlsode_source_file}
    TARGET_DIRECTORY ${PROJECT_NAME}
    PROPERTIES COMPILE_OPTIONS "${fortran_legacy_compile_option}"
)
