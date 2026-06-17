add_subdirectory(plugins)
target_sources(${PROJECT_NAME} PRIVATE
    system.f90
    stepper.f90
    driver.f90
)
