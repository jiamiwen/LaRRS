cmake_minimum_required(VERSION 3.23)
project(LaRRS C ASM)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(LaRRS main.c)

set(LARRS_SRCS sign.c packing.c polyvec.c poly.c ntt.c reduce.c)
set(LARRS_FIPS202_SRCS ${LARRS_SRCS} symmetric-shake.c)
set(FIPS202_SRCS fips202.c)
set(TEST_LARRS_SRCS test/test_LaRRS.c randombytes.c)
set(TEST_VECTORS_SRCS test/test_vectors.c)
set(TEST_SPEED_SRCS test/test_speed.c test/speed_print.c test/cpucycles.c randombytes.c)


if(MSVC)
    add_compile_options(/nologo /O2 /W4 /wd4146 /wd4244)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
    add_compile_options(-Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith)
    add_compile_options(-O0 -fomit-frame-pointer)
    add_compile_options(-fno-stack-check)
endif()

add_library(fips202_ref ${FIPS202_SRCS})


# LaRRS 2
add_library(larrs2_ref ${LARRS_FIPS202_SRCS})
target_compile_definitions(larrs2_ref PUBLIC LARRS_MODE=2)
target_link_libraries(larrs2_ref INTERFACE fips202_ref)

add_executable(test_larrs2_ref ${TEST_LARRS_SRCS})
add_executable(test_vectors2_ref ${TEST_VECTORS_SRCS})
add_executable(test_speed2_ref ${TEST_SPEED_SRCS})

target_link_libraries(test_larrs2_ref larrs2_ref)
target_link_libraries(test_vectors2_ref larrs2_ref)
target_link_libraries(test_speed2_ref larrs2_ref)