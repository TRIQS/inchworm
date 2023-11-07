
add_executable(integral integral.cpp)
target_link_libraries(integral xfac::xfac)
target_link_libraries(integral ${PROJECT_NAME}_c)

add_executable(integral_single_element integral_single_element.cpp)
target_link_libraries(integral_single_element xfac::xfac)
target_link_libraries(integral_single_element ${PROJECT_NAME}_c)

add_executable(integral_single_element_so integral_single_element_so.cpp)
target_link_libraries(integral_single_element_so xfac::xfac)
target_link_libraries(integral_single_element_so ${PROJECT_NAME}_c)

add_executable(integral_single_element_one_extra integral_single_element_one_extra.cpp)
target_link_libraries(integral_single_element_one_extra xfac::xfac)
target_link_libraries(integral_single_element_one_extra ${PROJECT_NAME}_c)

add_executable(integral_single_element_vertex integral_single_element_vertex.cpp)
target_link_libraries(integral_single_element_vertex xfac::xfac)
target_link_libraries(integral_single_element_vertex ${PROJECT_NAME}_c)

add_executable(integral_single_element_full integral_single_element_full.cpp)
target_link_libraries(integral_single_element_full xfac::xfac)
target_link_libraries(integral_single_element_full ${PROJECT_NAME}_c)

add_executable(test_decompose_so test_decompose_so.cpp)
target_link_libraries(test_decompose_so xfac::xfac)
target_link_libraries(test_decompose_so ${PROJECT_NAME}_c)

add_executable(integral_naive integral_naive.cpp)
target_link_libraries(integral_naive xfac::xfac)
target_link_libraries(integral_naive ${PROJECT_NAME}_c)
