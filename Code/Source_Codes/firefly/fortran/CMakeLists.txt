# construct the version file
execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION_DESCRIPTION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

string(REGEX REPLACE "^v" "" GIT_VERSION_DESCRIPTION "${GIT_VERSION_DESCRIPTION}")

configure_file(
    ${PROJECT_SOURCE_DIR}/config/version.f90.in
    ${CMAKE_CURRENT_BINARY_DIR}/version.f90
    @ONLY
)



# build the library
add_library(${LIBRARY_NAME} SHARED)
target_sources(${LIBRARY_NAME} PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/version.f90
    geometry.f90
    tasks.f90
)
set_source_files_properties(tasks.f90 PROPERTIES Fortran_PREPROCESS ON)


set(firefly_modules_directory ${CMAKE_BINARY_DIR}/modules)
set_target_properties(${LIBRARY_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_property(TARGET ${LIBRARY_NAME} PROPERTY Fortran_MODULE_DIRECTORY ${firefly_modules_directory})

target_compile_options(${LIBRARY_NAME} PRIVATE "-O2")
target_include_directories(${LIBRARY_NAME} PUBLIC ${firefly_modules_directory})
target_link_libraries(${LIBRARY_NAME} PUBLIC FLARE::FLARE)


install(TARGETS ${LIBRARY_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
