#include <iostream>
#include <complex>
#include "mfem.hpp"

#include "../headers/customcurl.hpp"
#include "solver.hpp"

using namespace mfem;

std::complex<real_t> getMu(std::complex<real_t> mu, Ferrite &ferrite, real_t &f,Mesh *mesh, real_t Bpeak, std::complex<real_t> mu0);


int main(int argc, char *argv[]){
    const char *path = "../../mesh/square.msh";             // Path to the mesh

    real_t f = 500e3;



    Mesh *mesh = new Mesh(path, 1, 1);
    mesh->UniformRefinement();

    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300, 1.8e6);
    Ferrite ferrite = N30;

    const std::complex<real_t> j(0., 1.);  // imaginary number such that j² = -1 
    const std::complex<real_t> mu_0 = ferrite.mu / ((real_t)1. + 1/ferrite.fc * j * 2.*M_PI * f);

    std::complex<real_t> mu = mu_0;
    for (int i = 0; i< 1000; i++){
        std::complex new_mu = getMu(mu, ferrite, f, mesh, 10e-3, mu_0);
        // std::cout << "Mu_r = " << new_mu.real() << ", Mu_i = " << new_mu.imag() << std::endl;
        std::cout << std::abs(new_mu) << std::endl;
        mu = new_mu;
    }

    return 0;
}

std::complex<real_t> getMu(std::complex<real_t> mu, Ferrite &ferrite, real_t &f,Mesh *mesh, real_t Bpeak, std::complex<real_t> mu0) {
    std::complex<real_t> new_mu(0,0);
// ---------------- Material Properties --------------
    real_t rho = ferrite.rho;
    real_t sigma = ferrite.sigma;
    real_t eps = ferrite.eps;


    real_t height = 7.59e-3;
    real_t w = 5.3e-3;

    // --------------- Complex values coefficients --------
    real_t omega = 2*M_PI*f;
    const std::complex<real_t> j(0., 1.);  // imaginary number such that j² = -1 

    const std::complex<real_t> rho_eq = rho + (real_t)1. /(sigma + j * omega * eps); 
    const std::complex<real_t> j_omega_mu = j*omega*mu;

    std::cout << "|rho| = " << std::abs(rho_eq) << std::endl;
    // The value of the source NI is initialized to 1. As the problem is linear, we can later rescale the solution.
    // We are doing such a thing because the actual source is the magnetic flux, but we don't know the corresponding current
    real_t NI = 1;

    // Boundary conditions function
    auto bdr_func = [NI](const Vector &x)
    {
        real_t r = x[0];
        return NI /(2 * M_PI * r);
    };

    // Function used for the integration because of the cylindrical coordinates
    auto r_coeff_func = [](const Vector &x){
        return x[0];  
    };

    // Function for the "correcting" factor rho/r² 
    auto inv_r_square_func = [](const Vector &x){
        return (real_t)1./pow(x[0],2);
    };


    int order = 1;                      // Order of shape functions
    int dim = mesh->Dimension();        // Mesh dimension (2D)

    // Working with scalar values as it is an axisymetric problem
    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    // Since we are in the harmonic steady-state regime, we work with complex values
    ComplexGridFunction h(fespace);

    // Definition of the Boundary conditions  :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 1; // All borders have Dirichlet boundary conditions (dbc), thus all border attributes are set to 1
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // Getting the list of degrees of freedom with dbc

    // FunctionCoefficient for the boundary conditions
    FunctionCoefficient bdr_coeff(bdr_func);  // Boundary condition of NI/(2*pi*r) on the border
    real_t zero = 0.;
    ConstantCoefficient zero_coeff(zero);     // Used for the complexe value

    // Projecting the boundary conditions on the border elements of the GridFunction h
    h.ProjectBdrCoefficient(bdr_coeff, zero_coeff, dir_bdr); // NI/(2*pi*r) for real part and 0 for imaginary part
    
    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Hermitian convention for the hermitian product
    b.Assemble();  // LinearForm equals to 0 as there is no source term

    // Coefficients for the definition of the bilinear form

    FunctionCoefficient r_coeff(r_coeff_func);   // Coefficient r for the integral in cylindrical coordinates

    // rho_eq Coefficient separated in real and imaginary part 
    ConstantCoefficient rho_real(rho_eq.real()); // for the -div(rho grad(h)) terme
    ConstantCoefficient rho_imag(rho_eq.imag());
    
    // **** Coefficient for the integration in cylindrical coordinates 
    ProductCoefficient r_rho_real(rho_real, r_coeff);   // r * rho_r
    ProductCoefficient r_rho_imag(rho_imag, r_coeff);   // r * rho_i
    // *****

    // j*omega*mu coefficient separated in real and imaginary part
    ConstantCoefficient j_omega_mu_real(j_omega_mu.real());  // jw * mu_r
    ConstantCoefficient j_omega_mu_imag(j_omega_mu.imag());  // jw * mu_i
    
    // **** Coefficient j*omega*mu_eq*r for the integration in cylindrical coordinates 
    ProductCoefficient r_j_omega_mu_real(j_omega_mu_real, r_coeff);  // r * jw * mu_r
    ProductCoefficient r_j_omega_mu_imag(j_omega_mu_imag, r_coeff);  // r * jw * mu_i
    // ****

    // rho_eq/r² term for correcting the -div(rho grad(h)) into a rot(rot(h)) . e_theta
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);                  // 1/r
    ProductCoefficient rho_inv_square_coeff_r(inv_r_square_coeff, r_rho_real);  // r * rho_r/r²
    ProductCoefficient rho_inv_square_coeff_i(inv_r_square_coeff, r_rho_imag);  // r * rho_i/r²
    

    // As we are using complex values, it is not a bilinear form but a sesquilinear form
    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN); // Hermitian convention
    a->AddDomainIntegrator(new DiffusionIntegrator(r_rho_real), new DiffusionIntegrator(r_rho_imag)); //  -div(rho_eq grad(H_theta))
    a->AddDomainIntegrator(new MassIntegrator(r_j_omega_mu_real), new MassIntegrator(r_j_omega_mu_imag));             //  jw mu_eq H

    // rho_eq/r² term for "correcting" the -div(rho grad(h)) into a rot(rot(h)) . e_theta
    a->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r), new MassIntegrator(rho_inv_square_coeff_i));  
    a->Assemble();
    
    // Defining a preconditioner using the real part
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new DiffusionIntegrator(r_rho_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(r_j_omega_mu_real));
    pcOp->AddDomainIntegrator(new MassIntegrator(rho_inv_square_coeff_r));
    pcOp->Assemble();

    OperatorHandle A;
    Vector B, H; 
    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B); 

    // Setting a preconditioner 
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

    pc_r = new GSSmoother(*PCOp.As<SparseMatrix>());  // Using a Gauss-Seidel smoother

    real_t s = -1.0;
    pc_i = new ScaledOperator(pc_r, s);  // Same preconditioner for omaginary part, only multiplied by -1 

    BDP.SetDiagonalBlock(0, pc_r);
    BDP.SetDiagonalBlock(1, pc_i);
    BDP.owns_blocks = 1;

    // --------- Solver ----------
    GMRESSolver gmres;
    gmres.SetPreconditioner(BDP);  // Applying the preconditioner
    gmres.SetOperator(*A.Ptr()); 
    gmres.SetRelTol(1e-12);
    gmres.SetMaxIter(2000);
    gmres.SetPrintLevel(0);
    gmres.Mult(B, H);   // Solve the system

    // 
    a->RecoverFEMSolution(H,b,h);  // Get system solution in h ComplexGridFunction
   
    GridFunction h_r = h.real();  // Getting real and imaginary part
    GridFunction h_i = h.imag();

    // GridFunction used for computing integral on the domain (for the magnetic flux)
    GridFunction ones(fespace);
    ones = 1; // All values set to 1

    ConstantCoefficient one_coeff(1.);

    // Coefficients for computing the integral of B (magnetic flux)
    GridFunctionCoefficient H_r_coeff(&h_r);
    GridFunctionCoefficient H_i_coeff(&h_i);


    // LinearForm representing the integral of H
    LinearForm lf_flux(fespace);
    lf_flux.AddDomainIntegrator(new DomainLFIntegrator(one_coeff));
    lf_flux.Assemble();
    real_t flux_r = lf_flux(h_r);
    real_t flux_i = lf_flux(h_i);

    std::complex<real_t>flux = flux_r + j * flux_i;  // Flux of H
    flux *= mu;               // Flux of B

    real_t flux_abs = std::abs(flux);

    GridFunction B_average(fespace);    // Constant GridFunction containing the average targeted value of B
    B_average = Bpeak/sqrt(2);          // Working with RMS values

    real_t targeted_flux = lf_flux(B_average);

    real_t coeff_correction;
    coeff_correction = targeted_flux / flux_abs;

    // Adjusting h values such that the magnetic flux is equal to the desired value
    h_r *= coeff_correction;
    h_i *= coeff_correction;

    // Computing the coresponding NI 
    NI *= coeff_correction; // RMS value
    flux *= coeff_correction;

    new_mu = 10e-3*w*height/std::abs(mu0)/flux;

    // std::cout << "Corrective factor : " << coeff_correction << std::endl;
    // std::cout << "NI (peak value) = " << NI*sqrt(2) << std::endl;

    // Free memory
    delete fec;
    delete fespace;
    delete a;
    delete pcOp;
    return new_mu;
}