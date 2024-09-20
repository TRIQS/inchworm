file(GLOB TCI_HYB_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(tci_hyb ${TCI_HYB_FILES})
target_link_libraries(tci_hyb PRIVATE xfac::xfac)
target_link_libraries(tci_hyb PRIVATE ${PROJECT_NAME}_c)
target_link_libraries(tci_hyb PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(tci_hyb PRIVATE nvtx_wrapper)
target_link_libraries(tci_hyb PRIVATE MPI::MPI_CXX)

target_include_directories(tci_hyb PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

