# Exemples d'application de la bibliothèque MFEM

Ce projet contient divers exemples d'utilisation de la bibliothèque MFEM pour résoudre des problèmes de méthodes des éléments finis.

## Introduction

Les exemples couvrent divers aspects de l'utilisation de MFEM, y compris la configuration de problèmes, l'utilisation de solveurs parallèles, et l'intégration avec d'autres bibliothèques comme Hypre et Metis.

## Prérequis

- [MFEM](https://github.com/mfem/mfem)
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


# MFEM sous Windows

Deux options sont disponibles pour compiler MFEM sur Windows :

1. Visual Studio + CMake 

    Approche classique via l'interface graphique ou en ligne de commande avec cmake.

2. : MSYS2 + MinGW (Recommandé)

    Cette méthode est plus simple dès lors que l'on souhaite utiliser les bibliothèques [Hypre](https://github.com/hypre-space/hypre) et [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) car elle se rapproche d'un environnement Unix.

Ci-dessous se trouve un guide permettant d'installer et de compiler MFEM avec [Hypre](https://github.com/hypre-space/hypre) et [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) grâce à MSYS2 

## Installation et configuration de MSYS2
Téléchargez MSYS2 via la [page officielle](https://www.msys2.org/) puis installez le. Une fois cela fait, il est necessaire d'installer les élements necessaires à la compilation de mfem. Heureusement ceux-ci sont accéssibles via les paquets.

Rmq : Il faut lancer MSYS2 MINGW64

   ```bash
    pacman -S mingw-w64-x86_64-toolchain
    pacman -S mingw-w64-x86_64-cmake
    pacman -S make
    pacman -S mingw-w64-x86_64-msmpi  
    pacman -S mingw-w64-x86_64-hypre
    pacman -S mingw-w64-x86_64-metis
    pacman -S git  
   ``` 


## Compiler MFEM :
Vous pouvez obtenir le code source de MFEM en le téléchargeant depuis la [page officielle](https://mfem.org/) ou en le clonant via Git :
   ```bash
    git clone https://github.com/mfem/mfem.git
    cd mfem
   ```
   Une fois dans le dossier de mfem, vous pouvez générer les fichiers de compilation avec cmake. 
   ```bash
    mkdir build
    cd build
    cmake -G "MSYS Makefiles" -DMFEM_USE_MPI=YES ..
   ```
   Le parmamètre -DMFEM_USE_MPI=YES permet de préciser à cmake que l'on souhaite compiler mfem en mode parallèle, c'est à dire en prenant hypre et metis en compte.

   Une fois les fichiers de compilation génerés, il ne reste plus qu'a compiler
   ```bash
    make -j $(nproc)
    make install  # Peut nécessiter les droits administrateur
   ```
