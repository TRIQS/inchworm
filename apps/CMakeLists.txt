
add_executable(integral integral.cpp)
target_link_libraries(integral xfac::xfac)
target_link_libraries(integral ${PROJECT_NAME}_c)

add_executable(test_decompose_so test_decompose_so.cpp)
target_link_libraries(test_decompose_so xfac::xfac)
target_link_libraries(test_decompose_so ${PROJECT_NAME}_c)

add_executable(integral_naive integral_naive.cpp)
target_link_libraries(integral_naive xfac::xfac)
target_link_libraries(integral_naive ${PROJECT_NAME}_c)
