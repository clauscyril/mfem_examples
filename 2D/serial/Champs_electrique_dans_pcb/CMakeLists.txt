cmake_minimum_required(VERSION 3.10)

# set(CMAKE_BUILD_TYPE Debug)
project(PCB)

set(CMAKE_CXX_STANDARD 11)

find_package(MFEM REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main mfem) 
