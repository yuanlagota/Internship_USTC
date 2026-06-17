# generate automatic f2py interfaces
set(f2py_autowrap_input_files
    _model.f90
    _analysis.f90
    _tasks.f90
    _mmesh.f90
)
set(f2py_autowrap_output_files "")
foreach(input_file IN LISTS f2py_autowrap_input_files)
    set(output_file ${CMAKE_CURRENT_BINARY_DIR}/${input_file})

    add_custom_command(
        OUTPUT ${output_file}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/flare/f2py/autowrap ${output_file} ${CMAKE_CURRENT_SOURCE_DIR}/../fortran
        COMMENT "Generating automatic f2py interfaces ${input_file}..."
        VERBATIM
    )

    list(APPEND f2py_autowrap_output_files ${output_file})
endforeach()



# source files for f2py interface
set(f2py_source_files
    version.f90
    kinds.f90
    control.f90
    utils.f90
    geqdsk.f90
    model.f90
    analysis.f90
    mmesh.f90
    tasks.f90
    cli.f90
)
list(TRANSFORM f2py_source_files PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/flare/f2py/)



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



# automatic documentation procedures
function(autodoc procedure_type filename procedure_name)
    set(prefix "")
    if(ARGV3)
        set(prefix ${ARGV3})
    endif()

    execute_process(
        COMMAND sed -n "/${procedure_type} ${prefix}${procedure_name}(/,\${p;/  use /q}" ${filename}
        COMMAND sed "1d"
        COMMAND sed "$d"
        COMMAND sed "s/^  !//"
        OUTPUT_VARIABLE doc
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${procedure_name}.txt ${doc})
endfunction()

function(autodoc_function)
    autodoc("function" ${ARGN})
endfunction()

function(autodoc_subroutine)
    autodoc("subroutine" ${ARGN})
endfunction()

# include subdirectories for documentation
add_subdirectory(flare/analysis)
add_subdirectory(flare/tasks)



# build f2py library interface
set(f2py_module_name f2py)
set(f2py_module_c ${f2py_module_name}module.c)
set(f2py_wrappers ${f2py_module_name}-f2pywrappers2.f90)

add_custom_target(${project_name}_genpyf ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c})

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

set(f2py_target ${project_name}_f2py)
Python3_add_library(${f2py_target} MODULE WITH_SOABI
    ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}
    ${CMAKE_CURRENT_BINARY_DIR}/${f2py_wrappers}
    ${F2PY_INCLUDE_DIR}/fortranobject.c
    ${f2py_source_files}
)

target_link_libraries(${f2py_target} PRIVATE Python3::NumPy PUBLIC MOOSE::MOOSE ${PROJECT_NAME})
add_dependencies(${f2py_target} ${project_name}_genpyf)
target_include_directories(${f2py_target} PRIVATE ${F2PY_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(${f2py_target} PROPERTIES OUTPUT_NAME ${f2py_module_name})
set_target_properties(${f2py_target} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)



# install the Python package
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${project_name}
    DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}
    PATTERN "*.pyc" EXCLUDE
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.dist-info" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
    PATTERN "f2py" EXCLUDE
)

install(TARGETS ${f2py_target} DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/${project_name})
