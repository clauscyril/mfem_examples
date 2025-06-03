#include <iostream>
#include <mat.h>
#include <unordered_map>
#include <filesystem>

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
    std::ofstream mesh_file("./meshs/"+std::string(name)+".msh");
    std::cout << "Écriture dans : " << std::filesystem::absolute("meshs/"+std::string(name)+".msh") << std::endl;

    if (!std::filesystem::exists("./meshs")) {
        std::cout << "creating meshs folder" << std::endl;
        std::filesystem::create_directory("./meshs");
        std::cout << "meshs folder created" << std::endl;
    }

    if (!mesh_file) {
        std::cerr << "Error : Impossible to open folder for writing.\n";
        return;
    }

    // Mesh format compatible with gmsh and mfem
    mesh_file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";


    // 1. Collecter les PhysicalNames
    std::set<std::pair<int, std::string>> physical_names; // {id, name}

    // Surface patches
    for (auto fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
        auto patch = c3t3.surface_patch_index(*fit);
        int patch_id_1 = patch.first;
        int patch_id_2 = patch.second;

        patch_id_1 = patch_id_1 == 0 ? patch_id_2 : patch_id_1;
        patch_id_2 = patch_id_2 == 0 ? patch_id_1 : patch_id_2;

        if (patch_id_1 > 0)
            physical_names.insert({patch_id_1, "surface_" + std::to_string(patch_id_1)});
        if (patch_id_2 > 0)
            physical_names.insert({patch_id_2, "surface_" + std::to_string(patch_id_2)});
    }

    // Volumes
    for (auto cit = c3t3.triangulation().finite_cells_begin(); cit != c3t3.triangulation().finite_cells_end(); ++cit) {
        if (!c3t3.is_in_complex(cit)) continue;
        int label = c3t3.subdomain_index(cit);
        physical_names.insert({label, "volume_" + std::to_string(label)});
    }

    // 2. Écrire la section $PhysicalNames
    mesh_file << "$PhysicalNames\n";
    mesh_file << physical_names.size() << "\n";
    for (const auto& [id, name] : physical_names) {
        int dim = name.find("surface") != std::string::npos ? 2 : 3;
        mesh_file << dim << " " << id << " \"" << name << "\"\n";
    }
    mesh_file << "$EndPhysicalNames\n";


    // Wrtiting the vertices 
    std::unordered_map<Tr::Vertex_handle, int> vertex_indices;   
    int index = 1;

    // First, all the nodes are listed with an index and its coordinates
    mesh_file << "$Nodes\n";

    // // On écrit le nombre de noeuds (Distance entre deux itérateurs) 
    // // opti du tout car on parcours les itérateurs deux fois (pour compter le total et pour les ajouter)
    // // À voir s'il n'y a pas moyen d'avoir le nb total d'élement directement
    mesh_file << std::distance(c3t3.triangulation().finite_vertices_begin(),
                               c3t3.triangulation().finite_vertices_end()) << std::endl;

    // Loop for listing all nodes
    for (auto vit = c3t3.triangulation().finite_vertices_begin();            // Looping on every vertices iterator
         vit != c3t3.triangulation().finite_vertices_end(); ++vit) {
        vertex_indices[vit] = index;  // Node index
        const auto& p = vit->point(); // Point object containing the nodes coordinates
        mesh_file << index++ << " " << p.x() << " " << p.y() << " " << p.z() << "\n"; // "Writing index + coordinates in msh format"
    }
    mesh_file << "$EndNodes\n";  // All nodes are listed

    // Vector for the triangles (border elements) and tetrahedras 
    std::vector<std::tuple<int, std::vector<int>, int>> triangles;
    std::vector<std::tuple<int, std::vector<int>, int>> tets;
    
    // Looping on the borders elements (triangles) by looping on facet iterators
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

        patch_id_1 = patch_id_1==0 ? patch_id_2 : patch_id_1; 
        patch_id_2 = patch_id_2==0 ? patch_id_1 : patch_id_2;

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
        mesh_file << elem_index++ << " 2 2 " << patch_id_1<< " " << patch_id_2 << " "
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



// void Savemsh(C3t3 &c3t3, const char* name) {
//     const std::string filename = "./meshs/" + std::string(name) + ".msh";
//     if (!std::filesystem::exists("./meshs")) {
//         std::filesystem::create_directory("./meshs");
//     }
//     std::ofstream mesh_file(filename);
//     if (!mesh_file) {
//         std::cerr << "Error: Cannot open file for writing.\n";
//         return;
//     }
//     std::cout << "Écriture dans : " << std::filesystem::absolute(filename) << std::endl;

//     // === Écriture de l'en-tête msh ===
//     mesh_file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

//     // === Préparation des structures ===
//     std::unordered_map<Tr::Vertex_handle, int> vertex_indices;
//     std::set<int> volume_ids; // pour écrire les PhysicalNames
//     std::map<std::pair<int, int>, int> surface_patch_to_id; // map patch pair -> surface ID
//     std::set<std::array<int, 3>> written_faces; // pour éviter les doublons
//     int surface_offset = 100; // décalage pour éviter chevauchement avec volume IDs
//     int next_surface_id = surface_offset;
//     int vertex_index = 1;

//     // === Écriture des nœuds ===
//     mesh_file << "$Nodes\n";
//     mesh_file << std::distance(c3t3.triangulation().finite_vertices_begin(),
//                                c3t3.triangulation().finite_vertices_end()) << "\n";

//     for (auto vit = c3t3.triangulation().finite_vertices_begin();
//          vit != c3t3.triangulation().finite_vertices_end(); ++vit) {
//         vertex_indices[vit] = vertex_index;
//         const auto& p = vit->point();
//         mesh_file << vertex_index++ << " " << p.x() << " " << p.y() << " " << p.z() << "\n";
//     }
//     mesh_file << "$EndNodes\n";

//     // === Collecte des triangles (surface) ===
//     std::vector<std::tuple<int, std::vector<int>>> triangle_elements;
//     for (auto fit = c3t3.facets_in_complex_begin();
//          fit != c3t3.facets_in_complex_end(); ++fit) {
//         auto cell = fit->first;
//         int i = fit->second;

//         int patch_id_1 = c3t3.surface_patch_index(*fit).first;
//         int patch_id_2 = c3t3.surface_patch_index(*fit).second;

//         // nettoyage des 0 (zone "inexistante")
//         std::pair<int, int> patch_pair = std::minmax(patch_id_1, patch_id_2);
//         if (patch_pair.first == 0) patch_pair.first = patch_pair.second;
//         if (patch_pair.first == 0) patch_pair.first = 999;

//         // ID unique pour chaque surface
//         if (surface_patch_to_id.count(patch_pair) == 0) {
//             surface_patch_to_id[patch_pair] = next_surface_id++;
//         }
//         int surface_id = surface_patch_to_id[patch_pair];

//         std::vector<int> nodes;
//         for (int j = 0; j < 4; ++j)
//             if (j != i)
//                 nodes.push_back(vertex_indices[cell->vertex(j)]);

//         std::array<int, 3> face_sorted = {nodes[0], nodes[1], nodes[2]};
//         std::sort(face_sorted.begin(), face_sorted.end());

//         if (written_faces.insert(face_sorted).second) {
//             triangle_elements.emplace_back(surface_id, nodes);
//         }
//     }

//     // === Collecte des tétraèdres ===
//     std::vector<std::tuple<int, std::vector<int>>> tetra_elements;
//     for (auto cit = c3t3.triangulation().finite_cells_begin();
//          cit != c3t3.triangulation().finite_cells_end(); ++cit) {
//         if (!c3t3.is_in_complex(cit)) continue;

//         int label = c3t3.subdomain_index(cit);
//         volume_ids.insert(label);

//         std::vector<int> nodes;
//         for (int i = 0; i < 4; ++i)
//             nodes.push_back(vertex_indices[cit->vertex(i)]);

//         tetra_elements.emplace_back(label, nodes);
//     }

//     std::cout << "test" << std::endl;
//     // === Écriture des PhysicalNames === (MUST BE AFTER filling surface_patch_to_id & volume_ids)
//     mesh_file << "$PhysicalNames\n";
//     mesh_file << surface_patch_to_id.size() + volume_ids.size() << "\n";

//     for (const auto& [pair, id] : surface_patch_to_id) {
//         mesh_file << "2 " << id << " \"surface_" << pair.first << "_" << pair.second << "\"\n";
//     }
//     for (int vol : volume_ids) {
//         mesh_file << "3 " << vol << " \"volume_" << vol << "\"\n";
//     }
//     mesh_file << "$EndPhysicalNames\n";


//     std::cout << "test 2" << std::endl;

//     // === Écriture des éléments ===
//     mesh_file << "$Elements\n";
//     mesh_file << triangle_elements.size() + tetra_elements.size() << "\n";

//     int elem_index = 1;

//     // Triangles (surface)
//     for (const auto& [phys_id, nodes] : triangle_elements) {
//         mesh_file << elem_index++ << " 2 2 " << phys_id << " " << phys_id << " "
//                   << nodes[0] << " " << nodes[1] << " " << nodes[2] << "\n";
//     }

//     // Tétraèdres (volumes)
//     for (const auto& [phys_id, nodes] : tetra_elements) {
//         mesh_file << elem_index++ << " 4 2 " << phys_id << " " << phys_id << " "
//                   << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3] << "\n";
//     }

//     mesh_file << "$EndElements\n";
//     mesh_file.close();
// }


// void Savemsh(C3t3 &c3t3, const char* name) {
//     std::ofstream mesh_file("./meshs/" + std::string(name) + ".msh");
//     std::cout << "Écriture dans : " << std::filesystem::absolute("meshs/" + std::string(name) + ".msh") << std::endl;

//     if (!std::filesystem::exists("./meshs")) {
//         std::cout << "creating meshs folder" << std::endl;
//         std::filesystem::create_directory("./meshs");
//         std::cout << "meshs folder created" << std::endl;
//     }

//     if (!mesh_file) {
//         std::cerr << "Error : Impossible to open folder for writing.\n";
//         return;
//     }

//     // Mapping zone → PhysicalName ID (avec décalage)
//     std::map<int, int> patch_to_physical;
//     int patch_offset = 100; // Décalage pour PhysicalNames surface

//     // Collecte des zones surfaces pour PhysicalNames
//     for (auto fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
//         auto patch = c3t3.surface_patch_index(*fit);
//         int p1 = patch.first;
//         int p2 = patch.second;

//         if (p1 != 0 && patch_to_physical.count(p1) == 0)
//             patch_to_physical[p1] = patch_offset + p1;
//         if (p2 != 0 && patch_to_physical.count(p2) == 0)
//             patch_to_physical[p2] = patch_offset + p2;
//     }

//     // Collecte des zones volumes pour PhysicalNames volumes (si besoin)
//     // Ici je ne crée pas PhysicalNames volume, mais tu peux ajouter si nécessaire

//     // --- Écriture entête msh ---
//     mesh_file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

//     // --- Écriture PhysicalNames ---
//     mesh_file << "$PhysicalNames\n";
//     mesh_file << patch_to_physical.size() << "\n";
//     for (const auto& [patch_id, phys_id] : patch_to_physical) {
//         mesh_file << "2 " << phys_id << " \"surface_" << patch_id << "\"\n";
//     }
//     mesh_file << "$EndPhysicalNames\n";

//     // --- Écriture noeuds ---
//     std::unordered_map<Tr::Vertex_handle, int> vertex_indices;
//     int index = 1;

//     mesh_file << "$Nodes\n";
//     mesh_file << std::distance(c3t3.triangulation().finite_vertices_begin(),
//                                c3t3.triangulation().finite_vertices_end()) << "\n";

//     for (auto vit = c3t3.triangulation().finite_vertices_begin();
//          vit != c3t3.triangulation().finite_vertices_end(); ++vit) {
//         vertex_indices[vit] = index;
//         const auto& p = vit->point();
//         mesh_file << index++ << " " << p.x() << " " << p.y() << " " << p.z() << "\n";
//     }
//     mesh_file << "$EndNodes\n";

//     // --- Écriture éléments ---
//     int elem_index = 1;

//     // Comptage total éléments surface + volume
//     size_t nb_facets = std::distance(c3t3.facets_in_complex_begin(), c3t3.facets_in_complex_end());
//     size_t nb_cells = 0;
//     for (auto cit = c3t3.triangulation().finite_cells_begin();
//          cit != c3t3.triangulation().finite_cells_end(); ++cit) {
//         if (c3t3.is_in_complex(cit)) ++nb_cells;
//     }

//     mesh_file << "$Elements\n";
//     mesh_file << (nb_facets + nb_cells) << "\n";

//     // Éléments surfaces (triangles)
//     for (auto fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
//         auto cell = fit->first;
//         int i = fit->second;

//         int p1 = c3t3.surface_patch_index(*fit).first;
//         int p2 = c3t3.surface_patch_index(*fit).second;

//         std::vector<int> tags;
//         if (p1 != 0) tags.push_back(patch_to_physical[p1]);
//         if (p2 != 0 && p2 != p1) tags.push_back(patch_to_physical[p2]);

//         std::vector<int> nodes;
//         for (int j = 0; j < 4; ++j)
//             if (j != i)
//                 nodes.push_back(vertex_indices[cell->vertex(j)]);

//         mesh_file << elem_index++ << " 2 " << tags.size();
//         for (int t : tags) mesh_file << " " << t;
//         for (int n : nodes) mesh_file << " " << n;
//         mesh_file << "\n";
//     }

//     // Éléments volumes (tétrahèdres)
//     for (auto cit = c3t3.triangulation().finite_cells_begin();
//          cit != c3t3.triangulation().finite_cells_end(); ++cit) {
//         if (!c3t3.is_in_complex(cit)) continue;

//         int label = c3t3.subdomain_index(cit);
//         // Ici on utilise label directement comme Physical tag pour volume,
//         // tu peux ajouter un offset si tu veux (ex: +0 ou +1000)
//         mesh_file << elem_index++ << " 4 1 " << label << " ";
//         for (int i = 0; i < 4; ++i) {
//             mesh_file << vertex_indices[cit->vertex(i)] << " ";
//         }
//         mesh_file << "\n";
//     }

//     mesh_file.close();
// }


int main(int argc, char* argv[])
{

    const char* path = "../hand.mat";
    const char* name = "data";

    // const char* path = "../../Alvar_v16.mat";
    // const char* name = "voxelData";



    CGAL::Image_3 img = getMatlabImage(path, name);

    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(img);

    // Mesh criteria
    Mesh_criteria criteria(params::facet_angle(30).facet_size(4).facet_distance(6).
                                    cell_radius_edge_ratio(4).cell_size(4));

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::lloyd().odt().perturb().exude());


    // Output
    std::ofstream medit_file("meshs/hand_all_2.mesh");
    // // c3t3.output_to_medit(medit_file, false);    // Rebind to false | (use CGAL::IO::output_to_medit which is deprecated (should use CGAL::IO::write_MEDIT))
    // CGAL::IO::write_MEDIT(medit_file, c3t3, params::all_cells(false).all_vertices(true));    // Should be equivalent to the line before but it's not exactly the case (have to check between CGAL::IO::write_MEDIT and CGAL::IO::output_to_medit (deprecated))
    // medit_file.close();

    Savemsh(c3t3, "hand_all");

    return 0;
}