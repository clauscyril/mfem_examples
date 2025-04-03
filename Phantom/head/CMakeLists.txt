cmake_minimum_required(VERSION 3.10)

if(WIN32)
    # On cherche l'emplacement de Matlab 
    execute_process(
        COMMAND reg query "HKLM\\SOFTWARE\\MathWorks\\MATLAB" /s
        OUTPUT_VARIABLE MATLAB_REGISTRY_OUTPUT
        ERROR_QUIET
    )
    string(REGEX MATCH "[A-Z]:\.*" MATLAB_PATH ${MATLAB_REGISTRY_OUTPUT})
    string(REGEX REPLACE "[\r\n]+" " " MATLAB_PATH "${MATLAB_PATH}")
    string(STRIP "${MATLAB_PATH}" MATLAB_PATH)
endif()


set(CMAKE_BUILD_TYPE Debug)

project(MFEMExample)

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

set(CMAKE_CXX_STANDARD 11)

# Si les bibliothèques ne sont pas dans les variables d'environnement, il faut préciser leur chemin d'accès
find_package(MFEM REQUIRED) 
find_package(CGAL REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Boost 1.70 REQUIRED COMPONENTS container system)

# Chemin vers les fichiers d'en-tête MATLAB
set(MATLAB_INCLUDE_DIR "${MATLAB_PATH}/extern/include")

# Chemin vers les bibliothèques MATLAB
link_directories("${MATLAB_PATH}/bin/win64")

add_compile_options(-Wall -Wextra -Wpedantic)



# Définir les fichiers sources
set(SOURCES_HEART heart.cpp)
set(SOURCES_READMATDATA readMatData.cpp)
set(SOURCES_GENERATEMESH generateMesh.cpp)

include_directories(${MATLAB_INCLUDE_DIR})

add_executable(heart ${SOURCES_HEART})
target_link_libraries(heart
    CGAL::CGAL
    CGAL::CGAL_Core  # For CGAL exception
    Boost::container # For boost::container errors
    Boost::system    # Often required by Boost components
    gmp              # GMP library
    mpfr             # MPFR library
)

add_executable(readMatData ${SOURCES_READMATDATA})
target_link_libraries(readMatData mfem CGAL::CGAL gmp mpfr libmat.dll libmx.dll) 

add_executable(generateMesh ${SOURCES_GENERATEMESH})
target_link_libraries(generateMesh mfem CGAL::CGAL gmp mpfr libmat.dll libmx.dll) 