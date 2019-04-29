cmake_minimum_required(VERSION 2.8)
project(nfl CXX ASM-ATT)

include_directories(include include/nfl include/nfl/arch include/nfl/opt include/nfl/opt/arch include/nfl/prng include/nfl/primes)

file(GLOB NFLLIB_ASMS lib/prng/*.s)
file(GLOB NFLLIB_SRCS lib/prng/*.cpp lib/params/*.cpp)
file(GLOB NFLLIB_HDRS_0 include/nfl.hpp)
file(GLOB NFLLIB_HDRS_1 include/nfl/*.hpp)
file(GLOB NFLLIB_HDRS_2 include/nfl/arch/*.hpp)
file(GLOB NFLLIB_HDRS_3 include/nfl/opt/*.hpp)
file(GLOB NFLLIB_HDRS_4 include/nfl/opt/arch/*.hpp)
file(GLOB NFLLIB_HDRS_5 include/nfl/primes/*.hpp)
file(GLOB NFLLIB_PRNG_HDRS include/nfl/prng/*.hpp include/nfl/prng/*.h)

foreach(f ${NFLLIB_ASMS})
  set_property(SOURCE ${f} PROPERTY LANGUAGE ASM-ATT)
endforeach(f)

# shared library
add_library(nfllib SHARED ${NFLLIB_ASMS} ${NFLLIB_SRCS})
target_link_libraries(nfllib ${GMP_LIBRARY} ${MPFR_LIBRARY})
install(TARGETS nfllib LIBRARY DESTINATION lib)

# static library
add_library(nfllib_static STATIC ${NFLLIB_ASMS} ${NFLLIB_SRCS})
target_link_libraries(nfllib_static ${GMP_LIBRARY} ${MPFR_LIBRARY})
install(TARGETS nfllib_static ARCHIVE DESTINATION lib)

install(FILES ${NFLLIB_HDRS_0} DESTINATION include)
install(FILES ${NFLLIB_HDRS_1} DESTINATION include/nfl)
install(FILES ${NFLLIB_HDRS_2} DESTINATION include/nfl/arch)
install(FILES ${NFLLIB_HDRS_3} DESTINATION include/nfl/opt)
install(FILES ${NFLLIB_HDRS_4} DESTINATION include/nfl/opt/arch)
install(FILES ${NFLLIB_HDRS_5} DESTINATION include/nfl/primes)
install(FILES ${NFLLIB_PRNG_HDRS} DESTINATION include/nfl/prng)
