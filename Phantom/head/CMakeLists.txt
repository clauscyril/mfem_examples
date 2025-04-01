cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Release)
project(MFEMExample)

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

set(CMAKE_CXX_STANDARD 11)

find_package(MFEM REQUIRED)
find_package(CGAL REQUIRED)

add_executable(main heart.cpp)
target_link_libraries(main mfem CGAL::CGAL gmp mpfr) 