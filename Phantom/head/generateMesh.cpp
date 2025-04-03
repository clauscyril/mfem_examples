#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
// #include <CGAL/Implicit_to_labeled_point_function_wrapper.h>

#include <iostream>
#include <mat.h>  // 📌 Header MATLAB officiel
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
using namespace std;

// Classe pour interpréter les voxels
struct VoxelData {
    std::vector<std::vector<std::vector<int>>> voxels;
    double voxel_size;
    int target_value; // Organe à extraire

    // Fonction implicite : retourne -1 si voxel appartient à l'organe, 1 sinon
    double operator()(const Kernel::Point_3& p) const {
        int x = static_cast<int>(p.x() / voxel_size);
        int y = static_cast<int>(p.y() / voxel_size);
        int z = static_cast<int>(p.z() / voxel_size);
        
        if (x >= 0 && x < voxels.size() && 
            y >= 0 && y < voxels[0].size() && 
            z >= 0 && z < voxels[0][0].size()) {
            return (voxels[x][y][z] == target_value) ? -1.0 : 1.0;
        }
        return 1.0; // Extérieur
    }
};

VoxelData loadMatFile(const string& filename, const string& variable_name);


int main() {

    string filename = "../tete.mat";
    string variable_name = "data";  // Nom de la variable MATLAB

    VoxelData voxel_data = loadMatFile(filename, variable_name);

    voxel_data.target_value = 11;


    // 📌 Définition de la surface implicite pour cet organe
    CGAL::Implicit_surface_3<Kernel, VoxelData> surface(voxel_data, Kernel::Sphere_3(Kernel::Point_3(10, 10, 10), 10));

    // 📌 Génération du maillage
    CGAL::Surface_mesh_default_criteria_3<Kernel> criteria(30, 0.1, 0.1);
    Mesh mesh;
    CGAL::make_surface_mesh(mesh, surface, criteria, CGAL::Manifold_tag());

    // 📌 Sauvegarde en STL
    std::ofstream output("organ_11.stl");
    output << mesh;
    output.close();

    return 0;
}




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