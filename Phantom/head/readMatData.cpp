#include <iostream>
#include <mat.h>  // 📌 Header MATLAB officiel
#include <vector>

using namespace std;

// 📌 Structure pour stocker les voxels
struct VoxelData {
    vector<vector<vector<int>>> voxels;
    double voxel_size;
};

// 📌 Fonction pour lire un fichier .mat avec `mat.h`
VoxelData loadMatFile(const string& filename, const string& variable_name) {
    // 📌 Ouvre le fichier .mat
    MATFile *pmat = matOpen(filename.c_str(), "r");
    if (!pmat) {
        cerr << "Erreur : Impossible d'ouvrir " << filename << endl;
        exit(1);
    }

    // 📌 Charge la variable
    mxArray *array = matGetVariable(pmat, variable_name.c_str());
    if (!array) {
        cerr << "Erreur : Variable '" << variable_name << "' non trouvée" << endl;
        matClose(pmat);
        exit(1);
    }

    // 📌 Récupération des dimensions
    size_t rows = mxGetDimensions(array)[0];
    size_t cols = mxGetDimensions(array)[1];
    size_t slices = mxGetDimensions(array)[2];

    cout << "Matrice 3D trouvée : " << rows << "x" << cols << "x" << slices << endl;

    // 📌 Copie des données
    VoxelData voxel_data;
    voxel_data.voxel_size = 1.0; // Modifier selon le contexte
    voxel_data.voxels.resize(rows, vector<vector<int>>(cols, vector<int>(slices, 0)));

    cout << rows << ", " << cols << ", " << slices << endl;
    int8_T* data = static_cast<int8_T*>(mxGetData(array));
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            for (size_t k = 0; k < slices; k++) {
                // cout << static_cast<int>(data[i + j * rows + k * rows * cols]) << " , " << k << endl;
                voxel_data.voxels[i][j][k] = static_cast<int>(data[i + j * rows + k * rows * cols]); 
            }
        }
    }

    // Nettoyage et fermeture
    mxDestroyArray(array);
    matClose(pmat);

    return voxel_data;
}

// 📌 Test
int main() {
    string filename = "../../Alvar_v16.mat";
    string variable_name = "voxelData";  // Nom de la variable MATLAB

    VoxelData voxel_data = loadMatFile(filename, variable_name);

    // 📌 Vérification d'un voxel
    cout << "Valeur du voxel (10,10,10) : " << voxel_data.voxels[1300][300][1600] << endl;

    return 0;
}
