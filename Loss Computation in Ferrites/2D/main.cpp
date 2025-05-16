#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "get_curl.hpp"

#define M_PI 3.14159265358979323846

using namespace mfem;

real_t bdr_func(const Vector &x);
real_t inv_r_square_func(const Vector &x);

// paramètres 
real_t Ri = 4.8e-3;

real_t d_g = 10e-6;
real_t d_b = 0.6e-9;

real_t sigma_g = 25;
real_t sigma_b = 0.67e-3;

real_t eps_b = 33;

real_t delta_g = d_g/(d_g + d_b);
// real_t delta_b = 1 - delta_g;
real_t delta_b = d_b/(d_b + d_g);


real_t rho = delta_g/sigma_g;
// real_t sigma = 1.;
// real_t eps = 1.;
real_t mu = 10e3 * 4e-7 * M_PI;

real_t tau = 1e0;  // Parmètre déterminé par un fit 

real_t omega = 2*M_PI*1e0;

// Dimensions Torre
real_t height = 7.6e-3;
real_t w = 5.3e-3;

// Paramètre Conditions aux limites
real_t N = 1;
real_t I = 0.1;

real_t Bpeak = 10e-3;
real_t mu_g = delta_g * mu;
real_t Rm = Ri + w/2;

// Coefficients complexes
// j imaginaire pur 
const std::complex<real_t> j(0., 1.);

const std::complex<real_t> rho_eq = rho + delta_b /(sigma_b + j * omega * eps_b);  // On est "obliger" de caster 1 en real_t (Au moins metre 1.)
const std::complex<real_t> mu_eq = j * mu * omega / ((real_t)1. + tau * j * omega);

// L'équation général est la suivante : 
//  rot(rho_eq rot(H)) + mu_eq jw H = 0
// Cependant, le problème présente une axisymétrie  => H = H_theta(r,z) . e_theta
// On peut donc se contenter d'un problème 2D 
// Malheureusment, lorsque H = H_theta(r,z) . e_theta, on n'a pas l'égalité suivante : 
// rot(rot H) . e_theta = - div(grad(H_theta)) mais on a rot(rot H) . e_theta = - div(grad(H_theta)) - H_theta/r²

// L'équation a résoudre est donc : 
//  - div(rho_eq grad(H_theta)) - H_theta/r² + jw mu_eq H_theta = 0


int main(){
    const char *path = "../../mesh/square.msh";
    Mesh *mesh = new Mesh(path, 1, 1);
    mesh->UniformRefinement();
    // mesh->UniformRefinement();
    mesh->PrintInfo(std::cout);

    int order = 1;
    int dim = 2;

    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    ComplexGridFunction h(fespace);

    // Conditions aux limites :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 1; // En 2D, tous les bords sont imposés
    
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // On précise dans le fespace quelles sont les "lignes" de Dirichlets
    FunctionCoefficient bdr_coeff(bdr_func);
    real_t zero = 0.;
    ConstantCoefficient zero_coeff(zero);
    h.ProjectBdrCoefficient(bdr_coeff, zero_coeff, dir_bdr);


    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Convention du produit hermitien
    b.Assemble();


    ConstantCoefficient rho_real(rho_eq.real()); 
    ConstantCoefficient rho_imag(rho_eq.imag()); 
    
    ConstantCoefficient mu_real(mu_eq.real());
    ConstantCoefficient mu_imag(mu_eq.imag());

    // Terme -1/r² qui apparait en 2D
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);

    std::cout << "mu_real = " << mu_eq.real() << ", mu imag = " << mu_eq.imag() << std::endl;
    std::cout << "rho_real = " << rho_eq.real() << ", rho imag = " << rho_eq.imag() << std::endl;
    std::cin.get();

    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN);
    a->AddDomainIntegrator(new DiffusionIntegrator(rho_real), new DiffusionIntegrator(rho_imag));
    a->AddDomainIntegrator(new MassIntegrator(mu_real), new MassIntegrator(mu_imag));
    a->AddDomainIntegrator(new MassIntegrator(inv_r_square_coeff), new MassIntegrator(inv_r_square_coeff));  // On peut utiliser le même coeff pour la partie réelle et la imaginaire
    a->Assemble();

    // Définition du preconditionneur (Uniquement la partie réel choisie ici)
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new DiffusionIntegrator(rho_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(mu_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(inv_r_square_coeff));
    pcOp->Assemble();


    OperatorHandle A, Pc;
    Vector B, H; 

    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B);
    
    // Utilisation du préconditionneur `pc0p` dans le solveur
    pcOp->FormSystemMatrix(ess_tdof_list, Pc);
    GSSmoother M((SparseMatrix&)(*Pc));  

    PCG(*A, M, B, H, 1, 10000, 1e-12, 0.f);

    Array<int> blockOffsets;
    blockOffsets.SetSize(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = A->Height() / 2;
    blockOffsets[2] = A->Height() / 2;
    blockOffsets.PartialSum();

    BlockDiagonalPreconditioner BDP(blockOffsets);

    Operator * pc_r = NULL;
    Operator * pc_i = NULL;

    OperatorHandle PCOp;
    
    pcOp->SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
    pcOp->FormSystemMatrix(ess_tdof_list, PCOp);

    pc_r = new DSmoother(*PCOp.As<SparseMatrix>());

    real_t s = -1.0;
    
    pc_i = new ScaledOperator(pc_r, s);


    BDP.SetDiagonalBlock(0, pc_r);
    BDP.SetDiagonalBlock(1, pc_i);

    BDP.owns_blocks = 1;


    GMRESSolver gmres;
    gmres.SetPreconditioner(BDP);
    gmres.SetOperator(*A.Ptr());
    gmres.SetRelTol(1e-12);
    gmres.SetMaxIter(100000);
    gmres.SetPrintLevel(1);
    gmres.Mult(B, H);

    a->RecoverFEMSolution(H,b,h);

    GridFunction h_r = h.real();
    GridFunction h_i = h.imag();
    GridFunction H_module(fespace);
    GridFunction H_phase(fespace);

    for (int i = 0; i<mesh->GetNV(); i++) {
        std::complex<real_t> h(h_r(i), h_i(i));
        H_module(i) = std::abs(h);
        H_phase(i) = std::arg(h) * 180 / M_PI;
    }

    FiniteElementCollection *fec_E = new ND_FECollection(order, dim);
    FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec_E);

    // FiniteElementSpace fespace_E(mesh, fec_E);

    CurlCustomCoefficient curl_r_coeff(&h_r);
    CurlCustomCoefficient curl_i_coeff(&h_i);
    GridFunction J_r(fespace_E);
    GridFunction J_i(fespace_E);
    J_r.ProjectCoefficient(curl_r_coeff);
    J_i.ProjectCoefficient(curl_i_coeff);

    PowerLossCoefficient Power_loss_coef(fespace_E, rho_eq, curl_r_coeff, curl_i_coeff);
    LinearForm lf(fespace);
    lf.AddDomainIntegrator(new DomainLFIntegrator(Power_loss_coef));
    lf.Assemble();

    GridFunction ones(fespace);
    ones = 1;
    real_t P_loss_tot = lf(ones);
    real_t P_loss_by_vol_mean = P_loss_tot/w/height;

    std::cout << "NI = " << sqrt(2)*Bpeak/mu_g *M_PI * Rm << std::endl;

    std::cout << "Ploss(W/m) : " << P_loss_tot << std::endl;
    std::cout << "Ploss(W/m^3) : " << P_loss_by_vol_mean << std::endl;

    // real_t max_J_r = 0;
    // real_t mas_j_i = 0;
    // for 


    bool visualisation = 1;

    if (visualisation){
        // Visualisation GLVis
        char vishost[] = "localhost";
        int  visport   = 19916;
        socketstream sol_sock_r(vishost, visport);
        socketstream sol_sock_i(vishost, visport);
        // socketstream sol_sock_mod(vishost, visport);
        // socketstream sol_sock_phase(vishost, visport);
        sol_sock_r.precision(8);
        sol_sock_i.precision(8);
        // sol_sock_mod.precision(8);
        // sol_sock_phase.precision(8);
        sol_sock_r << "solution\n" << *mesh << h_r
                    << "window_title 'Solution: Real Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        sol_sock_i << "solution\n" << *mesh << h_i
                    << "window_title 'Solution: Imaginary Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        // sol_sock_mod << "solution\n" << *mesh << H_module
        //             << "window_title 'Solution: Module'" 
        //             << "pause\n" << "keys c\n" << std::flush;
        // sol_sock_phase<< "solution\n" << *mesh << H_phase
        //             << "window_title 'Solution: Phase'" 
        //             << "pause\n" << "keys c\n" << std::flush;

        socketstream sol_sock_j_r(vishost, visport);
        socketstream sol_sock_j_i(vishost, visport);
        sol_sock_j_r.precision(8);
        sol_sock_j_i.precision(8);
        sol_sock_j_r << "solution\n" << *mesh << J_r
                    << "window_title 'Solution: J Real Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
        sol_sock_j_i << "solution\n" << *mesh << J_i
                    << "window_title 'Solution: J Imaginary Part'" 
                    << "pause\n" << "keys c\n" << std::flush;
    }
    // ParaViewDataCollection paraview_dc("torre_solution", mesh);
    // paraview_dc.SetPrefixPath("paraview_output"); // Optionnel : dossier de sortie
    // paraview_dc.RegisterField("H_real", &h_r);
    // paraview_dc.RegisterField("H_imag", &h_i);
    // paraview_dc.SetLevelsOfDetail(order);
    // paraview_dc.SetHighOrderOutput(true);
    // paraview_dc.SetCycle(0);
    // paraview_dc.SetTime(0.0);
    // paraview_dc.Save();

    
    return 0;
}

real_t bdr_func(const Vector &x)
{
    
    real_t r = Ri + x[0];
    // real_t module = N*I/(2 * M_PI * r);

    // return N*I/(2 * M_PI * r);
    return sqrt(2)*Bpeak/mu *M_PI * Rm /(2 * M_PI * r);

}

real_t inv_r_square_func(const Vector &x){

    return -1/pow(x[0] + Ri,2);
}
