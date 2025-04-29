# Exemples d'application de la bibliothèque MFEM

Ce projet contient divers exemples d'utilisation de la bibliothèque MFEM pour résoudre des problèmes de méthodes des éléments finis.

## Introduction

Les exemples couvrent divers aspects de l'utilisation de MFEM, y compris la configuration de problèmes, l'utilisation de solveurs parallèles, et l'intégration avec d'autres bibliothèques comme Hypre et Metis.

## Prérequis

- MFEM 
- CMake (version 3.10 ou supérieure)
- Hypre (facultatif, pour les solveurs parallèles)
- Metis (facultatif, pour le partitionnement de graphes)

## Utilisation

Pour compiler les exemples, suivez ces étapes :

1. Clonez ce dépôt :
   ```bash
   git clone https://github.com/clauscyril/mfem_examples.git
   cd mfem_examples
    ```
2. Dirigez vous vers le chemin de l'exemple que vous souhaitez compiler 
    ```bash
   cd 2D/paralell
    ```
3. Compilez avec CMake, pour cela, créez un dossier build
   ```bash
    mkdir build
    cd build
    ```
    et executer CMake et compiler avec make
   ```bash
    cmake ..
    make
    ```
liste des commandes : 

pacman -Syu

pacman -S mingw-w64-x86_64-git mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake mingw-w64-x86_64-make

pacman -S mingw-w64-x86_64-msmpi

pacman -S mingw-w64-x86_64-hypre

pacman -S mingw-w64-x86_64-metis

pacman -S mingw-w64-x86_64-mfem
