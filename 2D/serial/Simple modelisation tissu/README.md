Les oscillations apparaissent dans les 3 cas suivants : 

-   lorsqu'on utilise des fonctions de tests d'ordre 3 ou supérieurs
-   Lorsqu'on affine une fois le maillage via UniformRefinement et qu'on utilise des fonctions de tests d'ordre 2 ou supérieur
-   Lorsqu'on affine 2 fois ou plus le maillage via UniformRefinement

Cela s'explique par un problème de définitions de conditions aux limites. On a des conditions aux limites de type Neumann homogène : J.n = 0 sur le bord. Cependant, dans les trois cas énoncés ci-dessus, des nouveaux nœuds où la normal de J ne corresponds pas à la normal du bord sont ajoutés ! 
Dans ce cas, le problème n'est pas correctement défini, ce qui génère des erreurs
