file(GLOB MODE_FILES "mode/*.cpp")

add_executable(tci_hyb tci_hyb.cpp ${MODE_FILES})
target_link_libraries(tci_hyb xfac::xfac)
target_link_libraries(tci_hyb ${PROJECT_NAME}_c)
target_link_libraries(tci_hyb OpenMP::OpenMP_CXX)

