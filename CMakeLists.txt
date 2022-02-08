cmake_minimum_required(VERSION 3.10.2)
project (BiddingGA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_MATH_OVERFLOW_ERROR_POLICY=ignore_error")

include_directories(${PROJECT_SOURCE_DIR}/include)
 
set(SOURCES
  ${PROJECT_SOURCE_DIR}/src/auctions/first_price.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/second_price.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/all_pay.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/common_value_signal.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/common_value_signal_second.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/common_value_endpoints.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/common_value_endpoints2.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/common_value_signal_endpoints.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/first_price_reverse.cc
  ${PROJECT_SOURCE_DIR}/src/auctions/first_price_2d.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/vector_ops.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/roulette_zeroed.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/keep_best.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/keep_commonest.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/roulette.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/tournament.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/tournament_mixed.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/tournament_poisson.cc  
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/elitism_decorator.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/ranked_weighted.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/selector/ranked_exponential.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/binary/bit_mutator.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/binary/byte_array_genotype.cc
  ${PROJECT_SOURCE_DIR}/src/genericga/binary/single_point_crossover.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/distribution.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/interval.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/function_ops.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/distribution_ops.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/order_statistic_ops.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/bounds.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/combination_generation.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/scatter.cc 
  ${PROJECT_SOURCE_DIR}/src/numericaldists/grid_multi.cc
  ${PROJECT_SOURCE_DIR}/src/numericaldists/grid.cc
  ${PROJECT_SOURCE_DIR}/src/biddingga/initializers_1d.cc
  ${PROJECT_SOURCE_DIR}/src/biddingga/initializers_2d.cc
  ${PROJECT_SOURCE_DIR}/src/biddingga/helpers.cc
  )

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(biddingga ${PROJECT_SOURCE_DIR}/bidding_ga.cc ${SOURCES})
add_executable(biddingga_real ${PROJECT_SOURCE_DIR}/bidding_ga_real.cc ${SOURCES})
add_executable(common_value ${PROJECT_SOURCE_DIR}/common_value.cc ${SOURCES})
add_executable(common_value_full ${PROJECT_SOURCE_DIR}/common_value_full.cc ${SOURCES})
add_executable(biddingga_2d ${PROJECT_SOURCE_DIR}/bidding_ga_2d.cc ${SOURCES})
add_executable(ga_sample ${PROJECT_SOURCE_DIR}/bidding_ga_func_fit_sample.cc ${SOURCES})
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



# Taken from GoogleTests Manual

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Download and unpack googletest at configure time
#configure_file(CMakeLists.txt.in ${CMAKE_BINARY_DIR}/googletest-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#  RESULT_VARIABLE result
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
#endif()
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
#  RESULT_VARIABLE result
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#  message(FATAL_ERROR "Build step for googletest failed: ${result}")
#endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
#set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
#add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
#                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
#                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
#if (CMAKE_VERSION VERSION_LESS 2.8.11)
#  include_directories("${gtest_SOURCE_DIR}/include")
#endif()
set(TEST_SOURCES
  ${PROJECT_SOURCE_DIR}/tests_main.cc
  ${PROJECT_SOURCE_DIR}/test/auctions/first_price_tests.cc
  ${PROJECT_SOURCE_DIR}/test/numericaldists/distribution_ops_tests.cc
  ${PROJECT_SOURCE_DIR}/test/numericaldists/function_ops_tests.cc
  ${PROJECT_SOURCE_DIR}/test/numericaldists/order_statistic_ops_tests.cc
  ${PROJECT_SOURCE_DIR}/test/numericaldists/interval_tests.cc
  ${PROJECT_SOURCE_DIR}/test/numericaldists/combination_generation_tests.cc
  )

#add_executable(tests_main tests_main.cc ${SOURCES} ${TEST_SOURCES})
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include(GoogleTest)

set(TEST_MAIN ${PROJECT_SOURCE_DIR}/tests_main.cc)

add_executable(tests_main ${TEST_MAIN} ${SOURCES} ${TEST_SOURCES})
gtest_discover_tests(tests_main)


target_link_libraries(tests_main gtest gmock gtest_main pthread)
add_test(NAME    tests_main 
         COMMAND tests_main)
