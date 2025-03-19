cmake_minimum_required(VERSION 3.10)

# set(CMAKE_BUILD_TYPE Debug)
project(MFEMExample)

set(CMAKE_CXX_STANDARD 11)

find_package(MFEM REQUIRED)

add_executable(main main.cpp)
# # target_link_libraries(premierprojet mfem ${MPI_LIBRARIES} HYPRE) 
target_link_libraries(main mfem) 