#include <iostream>
#include "mat.h"   // Permet de lire les fichiers .mat (Données MATLAB)

using namespace std;

int main(){

    // Ouvrir le fichier .mat
    const char *file = "D:/Documents/projets/Alvar/Alvar_v16.mat";
    MATFile *pmat = matOpen(file, "r");
    const char *name;
    int i;

    
    if (pmat == nullptr) {
        cout << "Erreur lors de l'ouverture du fichier .mat" << endl;
        return 1;
    }

    // Obtenir la liste des variables dans le fichier
    int ndir = 0;
    char **dir = matGetDir(pmat, &ndir);
    if (dir == nullptr) {
        std::cerr << "Erreur lors de la récupération de la liste des variables" << std::endl;
        matClose(pmat);
        return 1;
    }
    // Afficher les noms des variables
    std::cout << "Liste des variables dans le fichier .mat :" << std::endl;
    for (i = 0; i < ndir; ++i) {
        std::cout << dir[i] << std::endl;
        // mxFree(dir[i]); // Libérer la mémoire allouée pour chaque nom de variable
    }
    mxFree(dir); // Libérer la mémoire allouée pour le tableau de noms de variables

    // Fermer le fichier .mat
    matClose(pmat);

    pmat = matOpen(file, "r");
    if (pmat == nullptr) {
    printf("Error reopening file %s\n", file);
    return(1);
    }


    mxArray *pa;

    /* Get headers of all variables */
    // printf("\nExamining the header for each variable:\n");
    for (i=0; i < ndir; i++) {
    pa = matGetNextVariableInfo(pmat, &name);
    if (pa == NULL) {
    // printf("Error reading in file %s\n", file);
    return(1);
    }
    /* Diagnose header pa */
    // printf("According to its header, array %s has %d dimensions\n",
    //     name, mxGetNumberOfDimensions(pa));
    if (mxIsFromGlobalWS(pa))
        printf("  and was a global variable when saved\n");
    else
        printf("  and was a local variable when saved\n");
    mxDestroyArray(pa);
    }

   /* Reopen file to read in actual arrays. */
   if (matClose(pmat) != 0) {
    // printf("Error closing file %s\n",file);
    return(1);
    }
    pmat = matOpen(file, "r");
    if (pmat == NULL) {
    // printf("Error reopening file %s\n", file);
    return(1);
    }

    /* Read in each array. */
    printf("\nReading in the actual array contents:\n");
    for (i=0; i<ndir; i++) {
        pa = matGetNextVariable(pmat, &name);
        if (pa == NULL) {
        printf("Error reading in file %s\n", file);
        return(1);
        } 
        /*
        * Diagnose array pa
        */
        // printf("According to its contents, array %s has %d dimensions\n",
        //     name, mxGetNumberOfDimensions(pa));
        if (mxIsFromGlobalWS(pa))
    printf("  and was a global variable when saved\n");
        else
    printf("  and was a local variable when saved\n");
        mxDestroyArray(pa);
    }

    return 0;
}

/*
 * MAT-file diagnose program
 *
 * See the MATLAB API Guide for compiling information.
 *
 * Calling syntax:
 *
 *   matdgns <matfile>
 *
 * It will diagnose the MAT-file named <matfile>.
 *
 * This program demonstrates the use of the following functions:
 *
 *  matClose
 *  matGetDir
 *  matGetNextVariable
 *  matGetNextVariableInfo
 *  matOpen
 *
 * Copyright 1984-2003 The MathWorks, Inc.
 */
// #include <stdio.h>
// #include <stdlib.h>
// #include "mat.h"

// int diagnose(const char *file) {
//     MATFile *pmat;
//     const char **dir;
//     const char *name;
//     int	  ndir;
//     int	  i;
//     mxArray *pa;

//     printf("Reading file %s...\n\n", file);

//     /*
//     * Open file to get directory
//     */
//     pmat = matOpen(file, "r");
//     if (pmat == NULL) {
//     printf("Error opening file %s\n", file);
//     return(1);
//     }

//     /*
//     * get directory of MAT-file
//     */
//     dir = (const char **)matGetDir(pmat, &ndir);
//     if (dir == NULL) {
//     printf("Error reading directory of file %s\n", file);
//     return(1);
//     } else {
//     printf("Directory of %s:\n", file);
//     for (i=0; i < ndir; i++)
//         printf("%s\n",dir[i]);
//     }
//     mxFree(dir);

//     /* In order to use matGetNextXXX correctly, reopen file to read in headers. */
//     if (matClose(pmat) != 0) {
//     printf("Error closing file %s\n",file);
//     return(1);
//     }
//     pmat = matOpen(file, "r");
//     if (pmat == NULL) {
//     printf("Error reopening file %s\n", file);
//     return(1);
//     }

//     /* Get headers of all variables */
//     printf("\nExamining the header for each variable:\n");
//     for (i=0; i < ndir; i++) {
//     pa = matGetNextVariableInfo(pmat, &name);
//     if (pa == NULL) {
//     printf("Error reading in file %s\n", file);
//     return(1);
//     }
//     /* Diagnose header pa */
//     printf("According to its header, array %s has %d dimensions\n",
//         name, mxGetNumberOfDimensions(pa));
//     if (mxIsFromGlobalWS(pa))
//         printf("  and was a global variable when saved\n");
//     else
//         printf("  and was a local variable when saved\n");
//     mxDestroyArray(pa);
//     }

//     /* Reopen file to read in actual arrays. */
//     if (matClose(pmat) != 0) {
//     printf("Error closing file %s\n",file);
//     return(1);
//     }
//     pmat = matOpen(file, "r");
//     if (pmat == NULL) {
//     printf("Error reopening file %s\n", file);
//     return(1);
//     }

//     /* Read in each array. */
//     printf("\nReading in the actual array contents:\n");
//     for (i=0; i<ndir; i++) {
//         pa = matGetNextVariable(pmat, &name);
//         if (pa == NULL) {
//         printf("Error reading in file %s\n", file);
//         return(1);
//         } 
//         /*
//         * Diagnose array pa
//         */
//         printf("According to its contents, array %s has %d dimensions\n",
//             name, mxGetNumberOfDimensions(pa));
//         if (mxIsFromGlobalWS(pa))
//     printf("  and was a global variable when saved\n");
//         else
//     printf("  and was a local variable when saved\n");
//         mxDestroyArray(pa);
//     }

//     if (matClose(pmat) != 0) {
//         printf("Error closing file %s\n",file);
//         return(1);
//     }
//     printf("Done\n");
//     return(0);
// }

// int main(int argc, char **argv)
// {

//     int result;
//     const char* path = "D:/Documents/projets/Alvar/Alvar_v16.mat";
//     if (argc > 1)
//     result = diagnose(path);
//     else{
//     result = 0;
//     printf("Usage: matdgns <matfile>");
//     printf(" where <matfile> is the name of the MAT-file");
//     printf(" to be diagnosed\n");
//     }

//     return (result==0)?EXIT_SUCCESS:EXIT_FAILURE;

// }

