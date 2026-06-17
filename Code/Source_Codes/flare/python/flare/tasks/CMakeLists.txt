set(fortran_tasks_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/../../../fortran/tasks)

file(GLOB fortran_tasks_source_files "${fortran_tasks_source_dir}/*.f90")
foreach(filename IN LISTS fortran_tasks_source_files)
    get_filename_component(basename ${filename} NAME_WE)
    autodoc_subroutine(${filename} ${basename} "flare_task_")
endforeach()

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    DESTINATION ${CMAKE_INSTALL_PYTHON_LIBDIR}/${project_name}/tasks
    FILES_MATCHING PATTERN *.txt
    PATTERN "CMakeFiles" EXCLUDE
)
