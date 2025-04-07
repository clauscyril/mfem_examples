#include <iostream>
#include <fstream>
#include <vector>

// #include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
// #include <CGAL/Implicit_surface_3.h>
// #include <CGAL/Surface_mesh_default_criteria_3.h>
// #include <CGAL/make_surface_mesh.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
 
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

// #include <CGAL/Mesh_3/config.h>
// #include <CGAL/Mesh_3/Implicit_surface_3.h>
// #include <CGAL/Mesh_3/make_surface_mesh.h>
// #include <CGAL/Mesh_3/Surface_mesher_generator.h>
// #include <CGAL/Mesh_3/Surface_mesher.h>

#include <mat.h>  // Header MATLAB officiel

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
    // 🔹 Charger le fichier .mat
    std::string filename = "../head.mat";
    std::string variable_name = "data";
    auto voxel_mat = loadMatFile(filename, variable_name);

    // 🔹 Préparer la structure CGAL
    VoxelData voxel_data;
    voxel_data.voxels = voxel_mat.voxels;
    voxel_data.voxel_size = 0.5116279069767442;  // mm (à ajuster selon ta donnée réelle)
    voxel_data.target_value = 11;  // Numéro d’organe à extraire

    // 🔹 Calcul de la bounding sphere automatiquement
    Kernel::Point_3 center(
        voxel_data.voxels.size() * voxel_data.voxel_size / 2.0,
        voxel_data.voxels[0].size() * voxel_data.voxel_size / 2.0,
        voxel_data.voxels[0][0].size() * voxel_data.voxel_size / 2.0
    );
    double radius = std::sqrt(3.0) * std::max({voxel_data.voxels.size(), voxel_data.voxels[0].size(), voxel_data.voxels[0][0].size()}) * voxel_data.voxel_size;

    // 🔹 Surface implicite pour CGAL
    CGAL::Implicit_surface_3<Kernel, VoxelData> surface(voxel_data, Kernel::Sphere_3(center, radius * radius));

    // Critères du maillage
    CGAL::Surface_mesh_default_criteria_3<Kernel> criteria(30, 0.5, 0.5); // Angle, taille max, distance approx

    //  Génération du mesh
    Mesh mesh;
    CGAL::make_surface_mesh(mesh, surface, criteria, CGAL::Manifold_tag());

    //  Sauvegarde STL
    std::ofstream output("organ_11.stl");
    output << mesh;
    output.close();

    std::cout << "Mesh exporté dans organ_11.stl" << std::endl;
    return 0;
}



//  Fonction pour lire un fichier .mat avec `mat.h`
VoxelData loadMatFile(const string& filename, const string& variable_name) {
    // Ouvre le fichier .mat
    MATFile *pmat = matOpen(filename.c_str(), "r");
    if (!pmat) {
        cerr << "Erreur : Impossible d'ouvrir " << filename << endl;
        exit(1);
    }

    // Charge la variable
    mxArray *array = matGetVariable(pmat, variable_name.c_str());
    if (!array) {
        cerr << "Erreur : Variable '" << variable_name << "' non trouvée" << endl;
        matClose(pmat);
        exit(1);
    }

    // Récupération des dimensions
    size_t rows = mxGetDimensions(array)[0];
    size_t cols = mxGetDimensions(array)[1];
    size_t slices = mxGetDimensions(array)[2];

    cout << "Matrice 3D trouvée : " << rows << "x" << cols << "x" << slices << endl;

    // Copie des données
    VoxelData voxel_data;
    voxel_data.voxel_size = 1.0; // Modifier selon le contexte
    voxel_data.voxels.resize(rows, vector<vector<int>>(cols, vector<int>(slices, 0)));

    cout << rows << ", " << cols << ", " << slices << endl;
    int8_T* data = static_cast<int8_T*>(mxGetData(array));    // Les données du voxel sont en int8 !
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