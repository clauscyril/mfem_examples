#include "solver.hpp"
#include <iostream>
#include <cmath>


using namespace mfem;

#define M_PI 3.14159265358979323846

// paramètres 
real_t Ri = 9.6e-3/2.0;

// // Paramètres Ferrite N30
real_t rho = 5.98e-2;
real_t sigma = 4.44e-1;
real_t eps = 2.48e-6;
real_t mu = 4300.0 * 4e-7 * M_PI;

// Paramètres Ferrite N87
// real_t rho = 4.24e-2;
// real_t sigma = 1.48e-1;
// real_t eps = 2.68e-6;
// real_t mu = 2200.0 *4e-7 * M_PI;


// Dimensions Torre
real_t height = 7.6e-3;
real_t w = 5.3e-3;

// // Paramètre Conditions aux limites
// real_t N = 1.;
// real_t I = 0.1;

real_t Bpeak = 10e-3;
real_t Rm = Ri + w/2;
real_t NI;

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

    Vector x; // Vector of the coordinates (r, z)
    T.Transform(ip, x); // Get the coordinates of the point ip in vector x, because ip has its coordinates in the element's referential
                        // Récupère les coordonnées du point ip dans le vector x (Car ip contient les coordonnées dans le référentiel d'un élément)

    real_t H_val = GridFunc->GetValue(T, ip);   

    Vector grad_H(2);
    GridFunc->GetGradient(T, grad_H); // On écrit les valeurs du gradient dans le vector grad_H.
                                      // Pas besoin de préciser le point ici car il s'agit du point attribué
                                      // dans l'ElementTransformation T (Voir définition de GetGradient)

    V[0] = -grad_H[1];                  // Caclul du rotationnel dans le plan (r,z) pour un vecteur de la forme :
    V[1] = grad_H[0] + H_val /(x[0]);   // H = H_theta(r,z) . e_theta   En ne prenant en entrée la GridFunction scalaire de H_theta
    // V[0] = -grad_H[1];                  // Caclul du rotationnel dans le plan (r,z) pour un vecteur de la forme :
    // V[1] = grad_H[0];   // H = H_theta(r,z) . e_theta   En ne prenant en entrée la GridFunction scalaire de H_theta
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
    return  (J_r_vect * J_r_vect + J_i_vect * J_i_vect) * rho_eq.real();
}

// Boundary conditions function
real_t bdr_func(const Vector &x)
{
    real_t r = x[0];
    // real_t r = Ri + w/2;
    // real_t module = N*I/(2 * M_PI * r);
    // return N*I/(2 * M_PI * r);
    // return sqrt(2)*Bpeak/mu * Rm /(2 * r);
    return NI /(2 * M_PI * r);
    // return sqrt(2)*Bpeak/mu * Rm /(2 * Rm);
}

// Function for the correcting factor rho/r²
real_t inv_r_square_func(const Vector &x){

    return 1/pow(x[0],2);
}

// Main function for computing the power loss
void GetPowerLossByFlux(const char* path, real_t fc, real_t fc_mu, real_t & P_loss_by_vol_mean, std::complex<real_t> &phi) {

    real_t omega = 2*M_PI*fc;
    real_t tau = 1./(2.*M_PI*fc_mu);  

    // Paramètres complexes 
    const std::complex<real_t> j(0., 1.);
    const std::complex<real_t> rho_eq = rho + (real_t)1. /(sigma + j * omega * eps); 
    const std::complex<real_t> mu_eq = mu / ((real_t)1. + tau * j * omega);

    const std::complex<real_t> j_omega_mu = j*omega*mu_eq;


    Mesh *mesh = new Mesh(path, 1, 1);
    int order = 1;                      // elements order : 1
    int dim = mesh->Dimension();        // Mesh dimension : 2

    // On travaille avec un scalaire, on définit donc un espace d'éléments H1
    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    // Since we are in the harmonic steady-state regime, we work with complex values
    ComplexGridFunction h(fespace);

    // //  Boundary conditions (no boundary condition) :
    Array<int> ess_tdof_list;
    ess_tdof_list = 0;

    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Hermitian convention for the hermitian product
    b.Assemble();  // LinearForm equals to 0

    // Coefficients for the definition of the bilinear form

    ConstantCoefficient rho_real(rho_eq.real()); // for the -div(rho grad(h)) terme
    ConstantCoefficient rho_imag(rho_eq.imag()); 

    ConstantCoefficient mu_real(mu_eq.real()); // for the -div(rho grad(h)) terme
    ConstantCoefficient mu_imag(mu_eq.imag()); 

    ConstantCoefficient j_omega_mu_real(j_omega_mu.real());  // for the jw mu h terme 
    ConstantCoefficient j_omega_mu_imag(j_omega_mu.imag());

    // rho_eq/r² terme for correcting the -div(rho grad(h)) into a rot(rot(h)) . e_theta
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);
    ProductCoefficient rho_inv_square_coeff_r(inv_r_square_coeff, rho_real);
    ProductCoefficient rho_inv_square_coeff_i(inv_r_square_coeff, rho_imag);

   
    // Définition de la forme bilinéaire complexe (forme sesquilinéaire)
    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN); // Convention du produit hermitien
    a->AddDomainIntegrator(new DiffusionIntegrator(rho_real), new DiffusionIntegrator(rho_imag)); //  -div(rho_eq grad(H_theta))
    a->AddDomainIntegrator(new MassIntegrator(j_omega_mu_real), new MassIntegrator(j_omega_mu_imag));             //  jw mu_eq H

    // rho_eq/r² terme for correcting the -div(rho grad(h)) into rot(rot(h)).e_theta
    a->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r), new MassIntegrator(rho_inv_square_coeff_i));  
    a->Assemble(); 
    

    ComplexLinearForm *c = new ComplexLinearForm(fespace, ComplexOperator::HERMITIAN);        // Constraint of the integral of the domain
    c->AddDomainIntegrator(new DomainLFIntegrator(mu_real), new DomainLFIntegrator(mu_imag));
    c->Assemble();

    OperatorHandle A, Pc;
    // SparseMatrix A;
    Vector B, H, Phi; 
    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B);  // Définition des matrices / vecteurs du problèmes
    
    // SparseMatrix *A_sparse = dynamic_cast<SparseMatrix *>(A.Ptr());
    // Operator *A2 = A.Ptr();
    // ConstrainedSolver S(*A2, A)
   

    // // 7. Extraction des résultats
    // Vector re_h(n/2), im_h(n/2);
    // for (int i = 0; i < n; i++)
    // {
    //     re_h(i) = sol(i);
    //     im_h(i) = sol(i + n);
    // }
    // std::cout << sol.Size() << ", "  << n << std::endl;
    // // std::complex<real_t> lambda(sol(n), sol(n + 1));
    // std::cout << "test" << std::endl;

    // std::cout << sol(n) << std::endl;

    // std::cout << "test" << std::endl;
    // // Reconstruction solution complexe
    // GridFunction h_re(fespace), h_im(fespace);
    // std::cout << "tretestsgv" << std::endl;
    // h_re.MakeRef(fespace, re_h, 0);
    // h_im.MakeRef(fespace, im_h, 0);


    // // Visualisation GLVis
    // std::cout << "test" << std::endl;
    // char vishost[] = "localhost";
    // int  visport   = 19916;
    // socketstream sol_sock_r(vishost, visport);
    // socketstream sol_sock_i(vishost, visport);
    // // socketstream sol_sock_mod(vishost, visport);
    // // socketstream sol_sock_phase(vishost, visport);
    // sol_sock_r.precision(8);
    // sol_sock_i.precision(8);
    // std::cout << "tes2t" << std::endl;
    // // sol_sock_mod.precision(8);
    // // sol_sock_phase.precision(8);
    // sol_sock_r << "solution\n" << *mesh << h_re
    //             << "window_title 'Solution: Real Part'" 
    //             << "pause\n" << "keys c\n" << std::flush;
    // sol_sock_i << "solution\n" << *mesh << h_im
    //             << "window_title 'Solution: Imaginary Part'" 
    //             << "pause\n" << "keys c\n" << std::flush;



    // system.SetBlock(3, 0, new SparseMatrix(ImC));
    // system.SetBlock(3, 1, new SparseMatrix(ReC));

    // // 5. Construire le second membre RHS
    // Vector rhs(offsets.Last());
    // rhs = 0.0;
    // for (int i = 0; i < n; i++)
    // {
    //     std::complex<double> bval = B(i);  // B = second membre complexe
    //     rhs(i) = bval.real();        // partie réelle
    //     rhs(i + n) = bval.imag();    // partie imaginaire
    // }
    // rhs(2 * n)     = Phi.real();    // contrainte réelle
    // rhs(2 * n + 1) = Phi.imag();    // contrainte imaginaire

    // // 6. Résolution
    // GMRESSolver gmres;
    // gmres.SetOperator(system);
    // gmres.SetRelTol(1e-12);
    // gmres.SetMaxIter(1000);
    // gmres.SetPrintLevel(1);
    // Vector sol(offsets.Last());
    // gmres.Mult(rhs, sol);

    // // 7. Extraction des résultats
    // Vector re_h(n), im_h(n);
    // for (int i = 0; i < n; i++)
    // {
    //     re_h(i) = sol(i);
    //     im_h(i) = sol(i + n);
    // }
    // std::complex<double> lambda(sol(2 * n), sol(2 * n + 1));

    // // Reconstruction solution complexe
    // GridFunction h_re(&fespace), h_im(&fespace);
    // h_re.MakeRef(fespace, re_h, 0);
    // h_im.MakeRef(fespace, im_h, 0);








    // // GridFunction utilisé pour calculer l'intégrale d'un coefficient sur le domaine
    // GridFunction ones(fespace);
    // ones = 1; // Toutes les valeurs sont définies comme étant égales à 1

    // // On récupère sous forme de coefficient les valeurs de H
    // GridFunctionCoefficient H_r_coeff(&h_r);
    // GridFunctionCoefficient H_i_coeff(&h_i);

    // GridFunction B_r(fespace);
    // GridFunction B_i(fespace);

    // B_r = h_r;
    // B_r *= mu_eq.real();

    // B_i = h_i;
    // B_i *= mu_eq.imag();

    // GridFunctionCoefficient B_r_coeff(&B_r);
    // GridFunctionCoefficient B_i_coeff(&B_i);

    // // On définit un linearForm pour intégrer H et donc pour obtenir le flux
    // LinearForm lf_r_flux(fespace);
    // lf_r_flux.AddDomainIntegrator(new DomainLFIntegrator(B_r_coeff));
    // lf_r_flux.Assemble();
    // real_t flux_r = lf_r_flux(ones);

    // // De même avec la partie imaginaire
    // LinearForm lf_i_flux(fespace);
    // lf_i_flux.AddDomainIntegrator(new DomainLFIntegrator(B_i_coeff));
    // lf_i_flux.Assemble();
    // real_t flux_i = lf_i_flux(ones);

    // flux = sqrt(flux_r * flux_r + flux_i * flux_i) / w / height;

    // // On affiche le flux obtenu
    // std::cout << "Flux magnétique = " << flux_r << " + " << flux_i << "j" << std::endl;
    // std::cout << "Flux magnétique (module) = " << flux << std::endl;

    // real_t coeff_correction = 7.5e-3/flux;
    // // real_t coeff_correction = 10e-3/flux;
    // h_r *= coeff_correction;
    // h_i *= coeff_correction;
    // // flux = sqrt(flux_r * flux_r + flux_i * flux_i) / w / height * coeff_correction;

    // // On définit un espace de Nedelec pour le courant J
    // FiniteElementCollection *fec_E = new ND_FECollection(order, dim);
    // FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec_E);

    // // Calcul du rotationnel de H, donc de J, sous la forme d'un coefficient
    // CurlCustomCoefficient curl_r_coeff(&h_r);
    // CurlCustomCoefficient curl_i_coeff(&h_i);
    
    // // Définition des GridFunction pour J
    // GridFunction J_r(fespace_E);
    // GridFunction J_i(fespace_E);

    // // On projette les valeurs des coefficients sur les GridFunction
    // J_r.ProjectCoefficient(curl_r_coeff);
    // J_i.ProjectCoefficient(curl_i_coeff);

    // // À partir des valeurs de J et de rho, on peut calcuer la densité de puissance 
    // PowerLossCoefficient Power_loss_coef(fespace_E, rho_eq, curl_r_coeff, curl_i_coeff);

    // // On intégre le densité de puissance sur le domaine
    // LinearForm lf(fespace);
    // lf.AddDomainIntegrator(new DomainLFIntegrator(Power_loss_coef));
    // lf.Assemble();
    // real_t P_loss_tot = lf(ones); // Puissance en W/m car on est en 2D

    // P_loss_by_vol_mean = P_loss_tot/w/height;  // Puissance moyenne en W/m^3

    // std::cout << "NI = " << sqrt(2)*Bpeak/mu *M_PI * Rm << std::endl;
    // std::cout << "Ploss(W/m) : " << P_loss_tot << std::endl;
    // std::cout << "Ploss(W/m^3) : " << P_loss_by_vol_mean << std::endl;

    // std::cout << "w = " << w << std::endl;
    // std::cout << "h = " << height << std::endl;
    // std::cout << "Re(rho) = " << rho_eq.real() << std::endl;

    // // On Projette la densité de puissance sur un GridFunction (utile seulement pour l'afficher)
    // GridFunction losses(fespace);
    // losses.ProjectCoefficient(Power_loss_coef);

    // bool visualisation = 0;  // 1 : plot activés
    //                          // 0 : plot desactivés

    // if (visualisation){
    //     // Visualisation GLVis
    //     char vishost[] = "localhost";
    //     int  visport   = 19916;
    //     socketstream sol_sock_r(vishost, visport);
    //     socketstream sol_sock_i(vishost, visport);
    //     // socketstream sol_sock_mod(vishost, visport);
    //     // socketstream sol_sock_phase(vishost, visport);
    //     sol_sock_r.precision(8);
    //     sol_sock_i.precision(8);
    //     // sol_sock_mod.precision(8);
    //     // sol_sock_phase.precision(8);
    //     sol_sock_r << "solution\n" << *mesh << h_r
    //                 << "window_title 'Solution: Real Part'" 
    //                 << "pause\n" << "keys c\n" << std::flush;
    //     sol_sock_i << "solution\n" << *mesh << h_i
    //                 << "window_title 'Solution: Imaginary Part'" 
    //                 << "pause\n" << "keys c\n" << std::flush;


    //     socketstream sol_sock_j_r(vishost, visport);
    //     socketstream sol_sock_j_i(vishost, visport);
    //     sol_sock_j_r.precision(8);
    //     sol_sock_j_i.precision(8);
    //     sol_sock_j_r << "solution\n" << *mesh << J_r
    //                 << "window_title 'Solution: J Real Part'" 
    //                 << "pause\n" << "keys c\n" << std::flush;
    //     sol_sock_j_i << "solution\n" << *mesh << J_i
    //                 << "window_title 'Solution: J Imaginary Part'" 
    //                 << "pause\n" << "keys c\n" << std::flush;
    //     // socketstream sol_sock_losses(vishost, visport);

    //     // sol_sock_losses.precision(8);

    //     // sol_sock_losses << "solution\n" << *mesh << losses
    //     //             << "window_title 'Solution: J Real Part'" 
    //     //             << "pause\n" << "keys c\n" << std::flush;
        
    // }
    // // ParaViewDataCollection paraview_dc("torre_solution", mesh);
    // // paraview_dc.SetPrefixPath("paraview_output"); // Optionnel : dossier de sortie
    // // paraview_dc.RegisterField("H_real", &h_r);
    // // paraview_dc.RegisterField("H_imag", &h_i);
    // // paraview_dc.SetLevelsOfDetail(order);
    // // paraview_dc.SetHighOrderOutput(true);
    // // paraview_dc.SetCycle(0);
    // // paraview_dc.SetTime(0.0);
    // // paraview_dc.Save();

    // // flux = sqrt(flux_r * flux_r + flux_i * flux_i) / w / height;

    // delete mesh;
    // delete fec;
    // delete fespace;
    // delete a;
    // delete fec_E;
    // delete fespace_E;
}


