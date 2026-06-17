# source files
set(f2py_source_files
    version.f90
    kinds.f90
    analysis.f90
    geometry.f90
    tools.f90
)
list(TRANSFORM f2py_source_files PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/moose/f2py/)



# set default installation directory
if(NOT CMAKE_INSTALL_PYTHON_LIBDIR)
    set(python_site_packages python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages)
    set(CMAKE_INSTALL_PYTHON_LIBDIR ${CMAKE_INSTALL_LIBDIR}/${python_site_packages} CACHE PATH "Installation directory for Python modules")
endif()



# grab the variables from a local Python installation
execute_process(
    COMMAND ${Python3_EXECUTABLE}
    -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)



# build f2py library interface
set(f2py_module_name f2py)
set(f2py_module_c ${f2py_module_name}module.c)
set(f2py_wrappers ${f2py_module_name}-f2pywrappers2.f90)

add_custom_target(genpyf ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c})

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c} ${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrappers}
    COMMAND ${Python3_EXECUTABLE}
            -m numpy.f2py
            ${f2py_source_files}
            -m ${f2py_module_name}
            --lower
    DEPENDS ${PROJECT_NAME} ${f2py_autowrap_output_files} ${f2py_source_files}
    COMMENT "Generating f2py wrappers"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

Python3_add_library(${f2py_module_name} MODULE WITH_SOABI
    ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}
    ${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrappers}
    ${F2PY_INCLUDE_DIR}/fortranobject.c
    ${f2py_source_files}
)

target_link_libraries(${f2py_module_name} PRIVATE Python3::NumPy PUBLIC ${PROJECT_NAME})
add_dependencies(${f2py_module_name} genpyf)
target_include_directories(${f2py_module_name} PRIVATE ${F2PY_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(${f2py_module_name} PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set_target_properties(${f2py_module_name} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)



# install the Python package
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${project_name}
    DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}
    PATTERN "*.pyc" EXCLUDE
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.dist-info" EXCLUDE
    PATTERN "f2py" EXCLUDE
)

install(TARGETS ${f2py_module_name} DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/${project_name})
