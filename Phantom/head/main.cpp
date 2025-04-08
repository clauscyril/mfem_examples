#include <iostream>
#include <mat.h>

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
        std::cerr << "Error : Impossible to open  " << filename << std::endl;
        exit(1);
    }

    // Exctract variable from file
    mxArray *array = matGetVariable(pmat, variableName);
    if (!array) {
        std::cerr << "Erreur : Variable '" << variableName << "' not found" << std::endl;
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
    _image *image = _createImage(rows, cols, slices, 1,     // ! Matlab order is not the same as CGAL
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


    for (size_t i=0; i<rows*cols*slices; i++) {
        if (set_brain.find(static_cast<int>(data[i])) != set_brain.end()) { // For now, only the brain is imported
            dst[i] = data[i];
        } else {
            dst[i] = 0;
        }
    }
    // _writeImage(image, "tete.inr");
    return CGAL::Image_3(image);
}

int main(int argc, char* argv[])
{

    const char* path = "../tete.mat";
    const char* name = "data";

    CGAL::Image_3 img = getMatlabImage(path, name);

    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(img);

    // Mesh criteria
    Mesh_criteria criteria(params::facet_angle(30).facet_size(10).facet_distance(2).
                                    cell_radius_edge_ratio(2).cell_size(3));

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::lloyd().odt().perturb().exude().manifold());
    // Output
    std::ofstream medit_file("out.mesh");
    CGAL::IO::write_MEDIT(medit_file, c3t3);
    medit_file.close();

    return 0;
}