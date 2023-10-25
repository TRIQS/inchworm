add_executable(integral integral.cpp)
target_link_libraries(integral xfac::xfac)
target_link_libraries(integral ${PROJECT_NAME}_c)