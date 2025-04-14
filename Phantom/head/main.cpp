#include <iostream>
#include <mat.h>
#include <unordered_map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
 
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
 
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/lloyd_optimize_mesh_3.h>
 
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
 
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
 
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
 
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
 
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
 
namespace params = CGAL::parameters;


CGAL::Image_3 getMatlabImage(const char *filename, const char *variableName){
    // Load matlab data from matlab file .mat
    MATFile *pmat = matOpen(filename, "r");
    if (!pmat) {
        std::cerr << "Error : Impossible to open : " << filename << std::endl;
        exit(1);
    }

    // Exctract variable from file
    mxArray *array = matGetVariable(pmat, variableName);
    if (!array) {
        std::cerr << "Erreur : Variable : '" << variableName << "' not found" << std::endl;
        matClose(pmat);
        exit(1);
    }

    // Convert data to uint8_t
    uint8_t* data = static_cast<uint8_t*>(mxGetData(array));

    // Get Voxel Size
    size_t rows = mxGetDimensions(array)[0];
    size_t cols = mxGetDimensions(array)[1];
    size_t slices = mxGetDimensions(array)[2];
    std::cout << rows << ", " << cols << ", " << slices << std::endl;
    // Initialize Image compatible with CGAL 
    _image *image = _createImage(rows, cols, slices, 1,     // 
                                1.f, 1.f, 1.f, 
                                1,                           // 1 = 8 bit, 2 = 16 bit (for more than 255 labels)
                                WK_FIXED, SGN_UNSIGNED);

    uint8_t *dst = (uint8_t *) image->data;
    
    std::set<uint8_t> set_brain;
    set_brain.insert(1);
    set_brain.insert(4);
    set_brain.insert(13);
    set_brain.insert(38);
    set_brain.insert(43);
    set_brain.insert(59);
    set_brain.insert(62);
    set_brain.insert(74);

    std::set<uint8_t> set_bones;
    set_bones.insert(11);
    // set_bones.insert(71);
    set_bones.insert(17);

    std::set<uint8_t> set_tooth;
    // set_tooth.insert(11);
    set_tooth.insert(71);

    std::set<uint8_t> set_tendons;
    set_tendons.insert(42);
    set_tendons.insert(51);
    set_tendons.insert(72);

    std::set<uint8_t> set_all;
    for (int i = 1;i<86; i++) {
        set_all.insert(i);
    }

    std::set<uint8_t> set_applied;

    set_applied = set_all;

    int count = 0;
    for (size_t i=0; i<rows*cols*slices; i++) {
        if (set_applied.find(static_cast<int>(data[i])) != set_applied.end()) { // For now, only the brain is imported
            dst[i] = data[i];
            count ++;
        } else {
            dst[i] = 0;
        }
    }

    std::cout << count << std::endl;
    // _writeImage(image, "brain.inr");
    return CGAL::Image_3(image);
}


void Savemsh(C3t3 &c3t3, const char* name){   // Have to check differences between string and char
    // Sauvegarde du mesh au format msh, compatible mfem
    std::ofstream mesh_file("meshs/"+std::string(name)+".msh");
    mesh_file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

    // Wrtiting the vertices 
    std::unordered_map<Tr::Vertex_handle, int> vertex_indices;   
    int index = 1;
    mesh_file << "$Nodes\n";

    // On écrit le nombre de noeuds (Distance entre deux itérateurs) 
    // pas opti du tout car on parcours les itérateurs deux fois (pour compter le total et pour les ajouter)
    // À voir s'il n'y a pas moyen d'avoir le nb total d'élement directement
    mesh_file << std::distance(c3t3.triangulation().finite_vertices_begin(),
                               c3t3.triangulation().finite_vertices_end()) << std::endl;

    
    for (auto vit = c3t3.triangulation().finite_vertices_begin();            // Looping on every vertices
         vit != c3t3.triangulation().finite_vertices_end(); ++vit) {

        vertex_indices[vit] = index;
        const auto& p = vit->point();
        mesh_file << index++ << " " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    mesh_file << "$EndNodes\n";

    

    std::vector<std::tuple<int, std::vector<int>, int>> triangles;
    std::vector<std::tuple<int, std::vector<int>, int>> tets;
    
    // on parcourt tous les triangles de bord
    for (auto fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
        auto cell = fit->first;
        int i = fit->second;
        int patch_id = c3t3.surface_patch_index(*fit).first;
    
        std::vector<int> nodes;
        for (int j = 0; j < 4; ++j)
            if (j != i) nodes.push_back(vertex_indices[cell->vertex(j)]);
        
        triangles.push_back({patch_id, nodes, 2}); // 2 = type Gmsh triangle
    }
    
    // On parcourt tous les tetraèdres
    for (auto cit = c3t3.triangulation().finite_cells_begin();
        cit != c3t3.triangulation().finite_cells_end(); ++cit) {
        if (!c3t3.is_in_complex(cit)) continue;
    
        int label = c3t3.subdomain_index(cit);
        std::vector<int> nodes;
        for (int i = 0; i < 4; ++i)
            nodes.push_back(vertex_indices[cit->vertex(i)]);
    
        tets.push_back({label, nodes, 4}); // 4 = type Gmsh tétraèdre
    }

    // 2. Write elements (surface + volume)
    int elem_index = 1;
    mesh_file << "$Elements\n";

    mesh_file << tets.size() + triangles.size() << std::endl;

    //  2.1 Surfaces (facets on boundary) 
    for (auto fit = c3t3.facets_in_complex_begin();
         fit != c3t3.facets_in_complex_end(); ++fit) {
        auto cell = fit->first;
        int i = fit->second;

        int patch_id_1 = c3t3.surface_patch_index(*fit).first; // label zone 1
        int patch_id_2 = c3t3.surface_patch_index(*fit).second; // label zone 2

        int patch_id = min(patch_id_1, patch_id_2); // Test (on ne prends pas les doublons)

        patch_id_1 = patch_id_1==0 ? 100 : patch_id_1; 
        patch_id_2 = patch_id_2==0 ? 100 : patch_id_2;

        patch_id = patch_id==0 ? 100 : patch_id;

        std::vector<int> nodes;

        for (int j = 0; j < 4; ++j) {
            if (j != i)
                nodes.push_back(vertex_indices[cell->vertex(j)]);
        }

        // Triangle: type 2
        // mesh_file << elem_index++ << " 2 2 " << patch_id_1 << " " << patch_id_1 << " "
        //     << nodes[0] << " " << nodes[1] << " " << nodes[2] << "\n";
        // mesh_file << elem_index++ << " 2 2 " << patch_id_2 << " " << patch_id_2 << " "
        //     << nodes[0] << " " << nodes[1] << " " << nodes[2] << "\n";
        mesh_file << elem_index++ << " 2 2 " << patch_id<< " " << patch_id << " "
            << nodes[0] << " " << nodes[1] << " " << nodes[2] << "\n";

    }


    // --- 2.2 Volumes (tetrahedra) ---
    for (auto cit = c3t3.triangulation().finite_cells_begin();
         cit != c3t3.triangulation().finite_cells_end(); ++cit) {
        if (!c3t3.is_in_complex(cit)) continue;

        int label = c3t3.subdomain_index(cit);
        mesh_file << elem_index++ << " 4 2 " << label << " " << label << " ";
        for (int i = 0; i < 4; ++i) {
            mesh_file << vertex_indices[cit->vertex(i)] << " ";
        }
        mesh_file << "\n";
    }


    mesh_file.close();
}


int main(int argc, char* argv[])
{

    const char* path = "../main.mat";
    const char* name = "data";

    // const char* path = "../../Alvar_v16.mat";
    // const char* name = "voxelData";



    CGAL::Image_3 img = getMatlabImage(path, name);

    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(img);

    // Mesh criteria
    Mesh_criteria criteria(params::facet_angle(30).facet_size(3).facet_distance(3).
                                    cell_radius_edge_ratio(2).cell_size(2));

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::lloyd().odt().perturb().exude());


    // Output
    // std::ofstream medit_file("meshs/hand_bones.mesh");
    // c3t3.output_to_medit(medit_file, false);    // Rebind to false | (use CGAL::IO::output_to_medit which is deprecated (should use CGAL::IO::write_MEDIT))
    // CGAL::IO::write_MEDIT(medit_file, c3t3, params::all_cells(false).all_vertices(true));    // Should be equivalent to the line before but it's not exactly the case (have to check between CGAL::IO::write_MEDIT and CGAL::IO::output_to_medit (deprecated))
    // medit_file.close();

    Savemsh(c3t3, "hand_all_2");

    return 0;
}