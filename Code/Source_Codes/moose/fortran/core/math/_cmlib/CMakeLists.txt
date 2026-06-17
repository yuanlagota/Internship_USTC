set(cmlib_source_files
    machcon.f90
    xerror.f90
    dbsplin.f90
    dtensbs.f90
    fftpkg.f90
    quadpkd.f90
)

target_sources(${PROJECT_NAME} PRIVATE ${cmlib_source_files})

set_source_files_properties(${cmlib_source_files}
    TARGET_DIRECTORY ${PROJECT_NAME}
    PROPERTIES COMPILE_OPTIONS "${fortran_legacy_compile_option}"
)
