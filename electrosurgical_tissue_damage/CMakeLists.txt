cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Release)

project(Elec_damage)

# set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

set(CMAKE_CXX_STANDARD 11)

# Si les bibliothèques ne sont pas dans les variables d'environnement, il faut préciser leur chemin d'accès
find_package(MFEM REQUIRED) 
# find_package(CGAL REQUIRED)
# find_package(CGAL REQUIRED COMPONENTS Core)
# find_package(Boost 1.70 REQUIRED COMPONENTS container system)

add_compile_options(-Wall -Wextra -Wpedantic)

add_executable(main main.cpp)
# target_link_libraries(main mfem CGAL::CGAL gmp mpfr) 
target_link_libraries(main mfem) 
