configure_file(inchworm-config.cmake.in inchworm-config.cmake @ONLY)
configure_file(inchworm-config-version.cmake.in inchworm-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/inchworm-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/inchworm-config-version.cmake
  DESTINATION lib/cmake/inchworm
)

install(EXPORT inchworm-targets NAMESPACE inchworm:: DESTINATION lib/cmake/inchworm)
