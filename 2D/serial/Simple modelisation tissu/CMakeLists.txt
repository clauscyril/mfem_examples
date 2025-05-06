cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Debug)
project(MFEMExample)

set(CMAKE_CXX_STANDARD 11)

find_package(MFEM REQUIRED)
# add_compile_options(-Wall -Wextra -Wpedantic)

# add_executable(main main.cpp)
# target_link_libraries(main mfem) 

# add_executable(cas_reel Cas_reel.cpp)
# target_link_libraries(cas_reel mfem) 

add_executable(cas_cond conduction_simple.cpp)
target_link_libraries(cas_cond mfem) 

# add_executable(test test.cpp)
# target_link_libraries(test mfem) 