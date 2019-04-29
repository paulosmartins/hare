cmake_minimum_required(VERSION 2.8)
project(hare CXX)

if(APPLE)
  cmake_policy(SET CMP0042 NEW)
endif()

#enable_testing()

include(CMakeCompilers.txt)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

add_subdirectory(nfl)

include_directories(include nfl/include)

file(GLOB HARE_HDRS_0 include/*.hpp)
file(GLOB HARE_HDRS_1 include/rns/*.hpp)
file(GLOB HARE_HDRS_2 include/hpr/*.hpp)

install(FILES ${HARE_HDRS_0} DESTINATION include)
install(FILES ${HARE_HDRS_1} DESTINATION include/rns)
install(FILES ${HARE_HDRS_2} DESTINATION include/hpr)

add_executable(main2 src/main2.cpp)
target_link_libraries(main2 nfllib_static ${GMP_LIBRARY} ${MPFR_LIBRARY})

