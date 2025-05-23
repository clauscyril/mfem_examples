#include "solver.hpp"
#include <iostream>
#include <cmath>

using namespace mfem;

#define M_PI 3.14159265358979323846

// paramètres 
real_t Ri = 9.6e-3/2.0;

// // Paramètres Ferrite N30
// real_t rho = 5.98e-2;
// real_t sigma = 4.44e-1;
// real_t eps = 2.48e-6;
// real_t mu = 4300.0 * 4e-7 * M_PI;

// // Paramètres Ferrite N87
real_t rho = 4.24e-2;
real_t sigma = 1.48e-1;
real_t eps = 2.68e-6;
real_t mu = 2200.0 *4e-7 * M_PI;


// Dimensions Torre
real_t height = 7.6e-3;
real_t w = 5.3e-3;

// // Paramètre Conditions aux limites
// real_t N = 1.;
// real_t I = 0.1;

real_t Bpeak = 10e-3;
real_t Rm = Ri + w/2;

// Constructeur du coefficient permmettant de calculer le rotationnel à partir d'un GridFunction
CurlCustomCoefficient::CurlCustomCoefficient (const GridFunction *gf)
    : VectorCoefficient((gf) ? gf -> FESpace() -> GetMesh() -> SpaceDimension() : 0)
{
   GridFunc = gf;
}

void CurlCustomCoefficient::SetGridFunction(const GridFunction *gf)
{
   GridFunc = gf; vdim = (gf) ?
                         gf -> FESpace() -> GetMesh() -> SpaceDimension() : 0;
}

void CurlCustomCoefficient::Eval(Vector &V, ElementTransformation &T,
                                           const IntegrationPoint &ip)
{

    Vector x; // Coordonnées
    T.Transform(ip, x); // Récupère les coordonnées du point ip dans le vector x (Car ip contient les coordonnées dans le référentiel d'un élément)

    real_t H_val = GridFunc->GetValue(T, ip);   

    Vector grad_H;
    GridFunc->GetGradient(T, grad_H); // On écrit les valeurs du gradient dans le vector grad_H.
                                      // Pas besoin de préciser le point ici car il s'agit du point attribué
                                      // dans l'ElementTransformation T (Voir définition de GetGradient)

    V[0] = -grad_H[1];                  // Caclul du rotationnel dans le plan (r,z) pour un vecteur de la forme :
    V[1] = grad_H[0] + H_val /(x[0]);   // H = H_theta(r,z) . e_theta   En ne prenant en entrée la GridFunction scalaire de H_theta
} 

PowerLossCoefficient::PowerLossCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> rho_eq_, CurlCustomCoefficient &J_r_, CurlCustomCoefficient &J_i_)
    : fespace(fespace_),
      J_r(J_r_), J_r_vect(fespace_->GetMesh()->SpaceDimension()),
      J_i(J_i_), J_i_vect(fespace_->GetMesh()->SpaceDimension()),
      rho_eq(rho_eq_) {}

real_t PowerLossCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
    J_r.Eval(J_r_vect, T, ip);
    J_i.Eval(J_i_vect, T, ip);
    return  (J_r_vect * J_r_vect + J_i_vect * J_i_vect) *rho_eq.real();
}

// Fonction qui calcule les valeurs des conditions aux limites
real_t bdr_func(const Vector &x)
{
    real_t r = x[0];
    // real_t r = Ri + w/2;
    // real_t module = N*I/(2 * M_PI * r);
    // return N*I/(2 * M_PI * r);
    // return sqrt(2)*Bpeak/mu * Rm /(2 * r);
    return sqrt(2)*Bpeak/mu * Rm /(2 * r);
}

real_t inv_r_square_func(const Vector &x){

    return 1/pow(x[0],2);
}

void GetPowerLoss(const char* path, real_t fc, real_t fc_mu, real_t & P_loss_by_vol_mean, real_t &flux) {

    real_t omega = 2*M_PI*fc;
    real_t tau = 1./(2.*M_PI*fc_mu);  // Parmètre déterminé par un fit 

    // Paramètres complexes 
    const std::complex<real_t> j(0., 1.);
    const std::complex<real_t> rho_eq = rho + (real_t)1. /(sigma + j * omega * eps);  // On est "obliger" de caster 1 en real_t (Au moins metre 1.)
    const std::complex<real_t> mu_eq = mu * j * omega / ((real_t)1. + tau * j * omega);
    // const std::complex<real_t> mu_eq = mu * j * omega;



    Mesh *mesh = new Mesh(path, 1, 1);
    mesh->UniformRefinement();
    mesh->UniformRefinement();
    // mesh->PrintInfo(std::cout);

    int order = 1;  // éléments d'ordre 1
    int dim = mesh->Dimension();

    // On travaille avec un scalaire, on définit donc un espace d'éléments H1
    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    // On est en régime harmonique, on travaille donc en valeurs complexes
    ComplexGridFunction h(fespace);

    // Conditions aux limites :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 1; // tous les bords sont imposés
    
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // On précise dans ess_tdof_list les éléments concernés par des conditions de Dirichlet
    
    // On définit les valeurs aux bords  
    FunctionCoefficient bdr_coeff(bdr_func);
    real_t zero = 0.;
    ConstantCoefficient zero_coeff(zero);
    // On projecte sur le GridFunction les valeurs imposées sur les éléments concernés
    h.ProjectBdrCoefficient(bdr_coeff, zero_coeff, dir_bdr); 
    // Forme linéaire complexe
    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Convention du produit hermitien
    b.Assemble();  // La forme linéaire est nulle 

    // Coefficients nécessaires pour la définitions de la forme bilinéaire
    ConstantCoefficient rho_real(rho_eq.real()); 
    ConstantCoefficient rho_imag(rho_eq.imag()); 

    ConstantCoefficient mu_real(mu_eq.real());
    ConstantCoefficient mu_imag(mu_eq.imag());

    // Terme rho_eq/r² qui apparrait en 2D
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);
    ProductCoefficient rho_inv_square_coeff_r(inv_r_square_coeff, rho_real);
    ProductCoefficient rho_inv_square_coeff_i(inv_r_square_coeff, rho_imag);

    // std::cout << "mu_real = " << mu_eq.real() << ", mu imag = " << mu_eq.imag() << std::endl;
    // std::cout << "rho_real = " << rho_eq.real() << ", rho imag = " << rho_eq.imag() << std::endl;

    // Définition de la forme bilinéaire complexe (forme sesquilinéaire)
    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN); // Convention du produit hermitien
    a->AddDomainIntegrator(new DiffusionIntegrator(rho_real), new DiffusionIntegrator(rho_imag)); // terme en -div(rho_eq grad(H_theta))
    a->AddDomainIntegrator(new MassIntegrator(mu_real), new MassIntegrator(mu_imag));             // terme en jw mu_eq H
    // terme supplémentaire rho_eq/r² H pour corriger la différence entre div(grad) et rot(rot()) en 2D
    a->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r), new MassIntegrator(rho_inv_square_coeff_i));  
    a->Assemble(); // Calcule les intégrales ajoutées
    
    // Définition du preconditionneur (Uniquement la partie réel choisie ici)
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new DiffusionIntegrator(rho_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(mu_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r));
    pcOp->Assemble();

    OperatorHandle A, Pc;
    Vector B, H; 
    std::cout << "ttestZ2" ;
    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B);  // Définition des matrices / vecteurs du problèmes
      
    // Utilisation du préconditionneur `pc0p` dans le solveur
    pcOp->FormSystemMatrix(ess_tdof_list, Pc);
    GSSmoother M((SparseMatrix&)(*Pc));  

    PCG(*A, M, B, H, 0, 10000, 1e-12, 0.f);

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


    // // Il y a plusieurs types de préconditionneurs disponibles :
    // pc_r = new DSmoother(*PCOp.As<SparseMatrix>());

    pc_r = new GSSmoother(*PCOp.As<SparseMatrix>());

    // // Si on utilise les calculs parallèles (ParMesh, etc...), on peut utiliser celui-ci :
    // HypreBoomerAMG *amg = new HypreBoomerAMG();
    // amg->SetOperator(*PCOp.As<HypreParMatrix>());
    // pc_r = amg;

    real_t s = -1.0;
    pc_i = new ScaledOperator(pc_r, s);  // On utilise le même preconditionneur pour la partie imaginaire, mais multiplié par -1

    BDP.SetDiagonalBlock(0, pc_r);
    BDP.SetDiagonalBlock(1, pc_i);

    BDP.owns_blocks = 1;

    GMRESSolver gmres;
    gmres.SetPreconditioner(BDP);
    gmres.SetOperator(*A.Ptr());
    gmres.SetRelTol(1e-12);
    gmres.SetMaxIter(10000);
    gmres.SetPrintLevel(0);
    gmres.Mult(B, H);

    a->RecoverFEMSolution(H,b,h);  // On récupère la solution du système dans le GridFunction h
    // On récupère la partie réélles et la partie imaginaire
    GridFunction h_r = h.real(); 
    GridFunction h_i = h.imag();

    // GridFunction utilisé pour calculer l'intégrale d'un coefficient sur le domaine
    GridFunction ones(fespace);
    ones = 1; // Toutes les valeurs sont définies comme étant égales à 1

    // On récupère sous forme de coefficient les valeurs de H
    GridFunctionCoefficient H_r_coeff(&h_r);
    GridFunctionCoefficient H_i_coeff(&h_i);

    GridFunction B_r(fespace);
    GridFunction B_i(fespace);

    B_r = h_r;
    B_r *= rho_eq.real();

    B_i = h_i;
    B_i *= rho_eq.imag();

    GridFunctionCoefficient B_r_coeff(&B_r);
    GridFunctionCoefficient B_i_coeff(&B_i);

    // On définit un linearForm pour intégrer H et donc pour obtenir le flux
    LinearForm lf_r_flux(fespace);
    lf_r_flux.AddDomainIntegrator(new DomainLFIntegrator(B_r_coeff));
    lf_r_flux.Assemble();
    real_t flux_r = lf_r_flux(ones);

    // De même avec la partie imaginaire
    LinearForm lf_i_flux(fespace);
    lf_i_flux.AddDomainIntegrator(new DomainLFIntegrator(B_i_coeff));
    lf_i_flux.Assemble();
    real_t flux_i = lf_i_flux(ones);

    // On affiche le flux obtenu
    std::cout << "Flux magnétique = " << flux_r << " + " << flux_i << "j" << std::endl;
    std::cout << "Flux magnétique (module) = " << sqrt(flux_r * flux_r + flux_i * flux_i) << std::endl;

    // On définit un espace de Nedelec pour le courant J
    FiniteElementCollection *fec_E = new ND_FECollection(order, dim);
    FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec_E);

    // Calcul du rotationnel de H, donc de J, sous la forme d'un coefficient
    CurlCustomCoefficient curl_r_coeff(&h_r);
    CurlCustomCoefficient curl_i_coeff(&h_i);
    
    // Définition des GridFunction pour J
    GridFunction J_r(fespace_E);
    GridFunction J_i(fespace_E);

    // On projette les valeurs des coefficients sur les GridFunction
    J_r.ProjectCoefficient(curl_r_coeff);
    J_i.ProjectCoefficient(curl_i_coeff);

    // À partir des valeurs de J et de rho, on peut calcuer la densité de puissance 
    PowerLossCoefficient Power_loss_coef(fespace_E, rho_eq, curl_r_coeff, curl_i_coeff);

    // On intégre le densité de puissance sur le domaine
    LinearForm lf(fespace);
    lf.AddDomainIntegrator(new DomainLFIntegrator(Power_loss_coef));
    lf.Assemble();
    real_t P_loss_tot = lf(ones); // Puissance en W/m car on est en 2D

    P_loss_by_vol_mean = P_loss_tot/w/height;  // Puissance moyenne en W/m^3

    std::cout << "NI = " << sqrt(2)*Bpeak/mu *M_PI * Rm << std::endl;
    std::cout << "Ploss(W/m) : " << P_loss_tot << std::endl;
    std::cout << "Ploss(W/m^3) : " << P_loss_by_vol_mean << std::endl;

    // On Projette la densité de puissance sur un GridFunction (utile seulement pour l'afficher)
    GridFunction losses(fespace);
    losses.ProjectCoefficient(Power_loss_coef);

    bool visualisation = 0;  // 1 : plot activés
                             // 0 : plot desactivés

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


        // socketstream sol_sock_j_r(vishost, visport);
        // socketstream sol_sock_j_i(vishost, visport);
        // sol_sock_j_r.precision(8);
        // sol_sock_j_i.precision(8);
        // sol_sock_j_r << "solution\n" << *mesh << J_r
        //             << "window_title 'Solution: J Real Part'" 
        //             << "pause\n" << "keys c\n" << std::flush;
        // sol_sock_j_i << "solution\n" << *mesh << J_i
        //             << "window_title 'Solution: J Imaginary Part'" 
        //             << "pause\n" << "keys c\n" << std::flush;
        // socketstream sol_sock_losses(vishost, visport);

        // sol_sock_losses.precision(8);

        // sol_sock_losses << "solution\n" << *mesh << losses
        //             << "window_title 'Solution: J Real Part'" 
        //             << "pause\n" << "keys c\n" << std::flush;
        
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

    flux = sqrt(flux_r * flux_r + flux_i * flux_i) / w / height;

}
