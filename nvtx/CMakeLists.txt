add_subdirectory(NVTX) 

file(GLOB sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(nvtx_wrapper ${sources})
target_include_directories(nvtx_wrapper PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(nvtx_wrapper PUBLIC nvtx3-cpp)
