cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Release)

project(MFEMExample)

# set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

set(CMAKE_CXX_STANDARD 17)

# Si les bibliothèques ne sont pas dans les variables d'environnement, il faut préciser leur chemin d'accès
find_package(CGAL REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Boost 1.70 REQUIRED COMPONENTS container system)


# find_library(GMP_LIBRARY NAMES gmp PATHS D:/MSYS2/mingw64/lib)
# find_library(MPFR_LIBRARY NAMES mpfr PATHS D:/MSYS2/mingw64/lib)

# target_link_libraries(MyTarget PRIVATE ${GMP_LIBRARY} ${MPFR_LIBRARY})

add_compile_options(-Wall -Wextra -Wpedantic)

# Définir les fichiers sources
set(SOURCES_HEART heart.cpp)

add_executable(heart ${SOURCES_HEART})
target_link_libraries(heart
    CGAL::CGAL
    CGAL::CGAL_Core  # For CGAL exception
    Boost::container # For boost::container errors
    Boost::system    # Often required by Boost components
    gmp              # GMP library
    mpfr             # MPFR library
    # ${GMP_LIBRARY}
    # ${MPFR_LIBRARY}
)

# target_include_directories(heart PRIVATE ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})