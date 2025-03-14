cmake_minimum_required(VERSION 3.10)

# set(CMAKE_BUILD_TYPE Debug)
project(MFEMExample)

set(CMAKE_CXX_STANDARD 11)

find_package(MFEM REQUIRED)


# find_package(MPI REQUIRED)

## Afficher les variables trouvées pour MPI
# message(STATUS "MPI_C_COMPILER: ${MPI_C_COMPILER}")
# message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
# message(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
# message(STATUS "MPI_LIBRARY: ${MPI_LIBRARY}")

# include_directories(${MPI_INCLUDE_PATH} /home/cyril/Documents/MFEM/hypre-2.32.0/src/hypre/include)

# # Ajouter le répertoire des bibliothèques de Hypre
# link_directories(${MPI_LIBRARY} /home/cyril/Documents/MFEM/hypre-2.32.0/src/hypre/lib)


add_executable(premierprojet /home/cyril/Documents/MFEM/PremierProjet/main.cpp)
# # target_link_libraries(premierprojet mfem ${MPI_LIBRARIES} HYPRE) 
target_link_libraries(premierprojet mfem) 

# add_executable(test_exc /home/cyril/Documents/MFEM/PremierProjet/test.cpp)
# target_link_libraries(test_exc mfem) 