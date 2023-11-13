file(GLOB MODE_FILES "mode/*.cpp")

add_executable(tci_inchworm_hyb tci_inchworm_hyb.cpp ${MODE_FILES})
target_link_libraries(tci_inchworm_hyb xfac::xfac)
target_link_libraries(tci_inchworm_hyb ${PROJECT_NAME}_c)

