#include <iostream>
#include <complex>
#include "mfem.hpp"

#define M_PI 3.14159265358979323846

using namespace mfem;

void bdr_func(const Vector &x, Vector &H);

// paramètres 
real_t Ri = 9.6e-3/2.0;

// Paramètres Ferrite N30
real_t rho = 5.98e-2;
real_t sigma = 4.44e-1;
real_t eps = 2.48e-6;
real_t mu = 4300.0 * 4e-7 * M_PI;

// // Paramètres Ferrite N87
// real_t rho = 4.24e-2;
// real_t sigma = 1.48e-1;
// real_t eps = 2.68e-6;
// real_t mu = 2200.0 *4e-7 * M_PI;



// Dimensions Torre
real_t height = 7.6e-3;
real_t w = 5.3e-3;

// // Paramètre Conditions aux limites
// real_t N = 1.;
// real_t I = 10;
real_t NI;

real_t Bpeak = 10e-3;
// real_t mu_g = delta_g * mu;
real_t Rm = Ri + w/2;





class PowerLossCoefficient : public mfem::Coefficient
{
private:
    const FiniteElementSpace *fespace;
    CurlGridFunctionCoefficient J_r;
    mutable mfem::Vector J_r_vect;
    CurlGridFunctionCoefficient J_i;
    mutable mfem::Vector J_i_vect;
    std::complex<real_t> rho_eq;

public:
    PowerLossCoefficient(const FiniteElementSpace *fespace_, std::complex<real_t> rho_eq_, CurlGridFunctionCoefficient &J_r_, CurlGridFunctionCoefficient &J_i_)
    : fespace(fespace_),
      J_r(J_r_), J_r_vect(fespace_->GetMesh()->SpaceDimension()),
      J_i(J_i_), J_i_vect(fespace_->GetMesh()->SpaceDimension()),
      rho_eq(rho_eq_) {}
    
    virtual real_t Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
    {
        Vector x;
        T.Transform(ip, x);
        J_r.Eval(J_r_vect, T, ip);
        J_i.Eval(J_i_vect, T, ip);
        real_t r = sqrt(x[0]*x[0] + x[1]*x[1]);
        return  (J_r_vect * J_r_vect + J_i_vect * J_i_vect) *rho_eq.real() * r;
    }

    
};

void source(const Vector &x, Vector &H) {
    real_t r = sqrt(x[0]*x[0] + x[1]*x[1]);
    real_t module = NI/ (2 * M_PI* r);
    H(0) = -x[1]/r * module;
    H(1) = x[0]/r * module;
    H(2) = 0.0;
};

    

real_t ComputePowerLoss(real_t fc, real_t fc_mu){

    real_t omega = 2*M_PI*fc;
    real_t tau = 1./(2.*M_PI*fc_mu);  // Parmètre déterminé par un fit 

    const std::complex<real_t> j(0., 1.);

    const std::complex<real_t> rho_eq = rho + (real_t)1. /(sigma + j * omega * eps);  // On est "obliger" de caster 1 en real_t (Au moins metre 1.)
    const std::complex<real_t> mu_eq = mu / ((real_t)1. + tau * j * omega);
    const std::complex<real_t> j_omega_mu_eq = mu_eq * j * omega;

    const char *path = "../../mesh/torre_carre.msh";
    Mesh *mesh = new Mesh(path, 1, 1);
    // mesh->PrintInfo(std::cout);

    int order = 1;
    int dim = mesh->Dimension();

    FiniteElementCollection *fec = new ND_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    ComplexGridFunction h(fespace);



    // Conditions aux limites :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    // dir_bdr[26] = 1; 
    dir_bdr = 1; 
    NI = sqrt(2)*Bpeak/std::abs(mu_eq) * M_PI * Rm;
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);  // On précise dans le fespace quelles sont les "lignes" de Dirichlets
    VectorFunctionCoefficient bdr_coeff(dim, bdr_func);
    Vector zero(dim);
    zero = 0.;
    VectorConstantCoefficient zero_coeff(zero);

    h.ProjectBdrCoefficientTangent(bdr_coeff, zero_coeff, dir_bdr);

    // VectorFunctionCoefficient source_coeff(dim, source);

    ComplexLinearForm b(fespace, ComplexOperator::HERMITIAN); // Convention du produit hermitien
    // b.AddDomainIntegrator(new VectorFEDomainLFIntegrator(source_coeff), NULL);
    b.Assemble();
    // b.Assemble();

    // ConstantCoefficient rho_real(sigma/(sigma*sigma + pow(eps*omega, 2))); 
    // ConstantCoefficient rho_imag(-omega/(sigma*sigma + pow(eps*omega, 2))); 
    ConstantCoefficient rho_real(rho_eq.real()); 
    ConstantCoefficient rho_imag(rho_eq.imag()); 
    

    // ConstantCoefficient mu_real(mu*omega*omega*tau/(1+pow(tau*omega,2)));
    // ConstantCoefficient mu_imag(mu*omega/(1+pow(tau*omega,2)));
    ConstantCoefficient mu_real(j_omega_mu_eq.real());
    ConstantCoefficient mu_imag(j_omega_mu_eq.imag());


    std::cout << "mu_real = " << mu_eq.real() << ", mu imag = " << mu_eq.imag() << std::endl;

    SesquilinearForm *a = new SesquilinearForm(fespace, ComplexOperator::HERMITIAN);
    a->AddDomainIntegrator(new CurlCurlIntegrator(rho_real), new CurlCurlIntegrator(rho_imag));
    a->AddDomainIntegrator(new VectorFEMassIntegrator(mu_real), new VectorFEMassIntegrator(mu_imag));
    a->Assemble();

    // Définition du preconditionneur (Uniquement la partie réel choisie ici)
    BilinearForm *pcOp = new BilinearForm(fespace);
    pcOp->AddDomainIntegrator(new CurlCurlIntegrator(rho_real));
    pcOp->AddDomainIntegrator(new VectorFEMassIntegrator(mu_real));
    pcOp->Assemble();


    OperatorHandle A, Pc;
    Vector B, H; 

    a->FormLinearSystem(ess_tdof_list, h, b, A, H, B);
    

    std::cout << "test" << std::endl;
    // Utilisation du préconditionneur `pc0p` dans le solveur
    pcOp->FormSystemMatrix(ess_tdof_list, Pc);
    GSSmoother M((SparseMatrix&)(*Pc));  


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
    gmres.SetMaxIter(10000);
    gmres.SetPrintLevel(0);
    gmres.Mult(B, H);

    a->RecoverFEMSolution(H,b,h);

    GridFunction h_r = h.real();
    GridFunction h_i = h.imag();

    // real_t Brms = 
    // real_t coeff_correction;
    // coeff_correction = Bpeak/sqrt(2)/B_rms;

        
    // // real_t coeff_correction = 10e-3/flux;
    // h_r *= coeff_correction;
    // h_i *= coeff_correction;

    CurlGridFunctionCoefficient J_r(&h_r);
    CurlGridFunctionCoefficient J_i(&h_i);

    GridFunction J_r_grid(fespace);
    GridFunction J_i_grid(fespace);

    J_r_grid.ProjectCoefficient(J_r);
    J_i_grid.ProjectCoefficient(J_i);

    std::cout << "test" << std::endl;

    PowerLossCoefficient Power_loss_coef(fespace, rho_eq, J_r, J_i);

    FiniteElementCollection *fec_H1 = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace_H1 = new FiniteElementSpace(mesh, fec_H1);


    LinearForm lf(fespace_H1);
    lf.AddDomainIntegrator(new DomainLFIntegrator(Power_loss_coef));
    lf.Assemble();
    std::cout << "test2" << std::endl;

    GridFunction ones(fespace_H1);
    ones = 1;
    real_t P_loss_tot = lf(ones);
    real_t vol = M_PI*height * (pow(Ri + w,2) - pow(Ri,2));
    real_t P_loss_by_vol_mean = P_loss_tot / vol;

    std::cout << "Ploss(W/m) : " << P_loss_tot << std::endl;
    std::cout << "Ploss(W/m^3) : " << P_loss_by_vol_mean << std::endl;
    

//    // Visualisation GLVis
//     char vishost[] = "localhost";
//     int  visport   = 19916;
//     socketstream sol_sock_i(vishost, visport);
//     socketstream sol_sock_grad(vishost, visport);
//     sol_sock_i.precision(8);
//     sol_sock_grad.precision(8);
//     sol_sock_i << "solution\n" << *mesh << J_r_grid
//                 << "window_title 'Solution: Real Part'" 
//                 << "pause\n" << "keys c\n" << std::flush;
//     sol_sock_grad << "solution\n" << *mesh << J_i_grid
//                 << "window_title 'Solution: Imaginary Part'" 
//                 << "pause\n" << "keys c\n" << std::flush;

    return P_loss_by_vol_mean;
}

void bdr_func(const Vector &x, Vector &H)
{
    // Cas Torre :
    real_t r = sqrt(x[0]*x[0] + x[1]*x[1]);
    real_t module =NI/(2 * M_PI*r);;
    H(0) = -x[1]/r * module;
    H(1) = x[0]/r * module;
    H(2) = 0.;


    // Cas Cube : 
    // real_t Ri = 1;
    // real_t r = Ri + x[0];
    // H(0) = 0;
    // H(1) = 0;
    // H(2) = N*I/(2 * M_PI * r);
}


int main() {
    // for (int j = 0; j < 5; j++){
    //     std::string name = "./data";
    //     name += std::to_string(j) + ".csv";
    //     std::cout << name << std::endl;
    //     std::ofstream data_file(name);
    //     data_file << "fc;Ploss\n";
    //     real_t fc_0 = 50e3;
    //     real_t fc = fc_0;
    //     real_t fc_end = 2e6;
    //     int N = 14;
    //     real_t fc_mu =  1e6;
    //     real_t fc_mu_end = 10e6;
    //     real_t delta_fc_mu = (fc_mu_end-fc_mu)/4;
    //     fc_mu = 1e6 + j*delta_fc_mu;
    //     // real_t delta_f = (fc_end - fc)/N;
    //     std::cout << fc_mu << std::endl;

    //     real_t u = log(fc_0);
    //     real_t u_end = log(fc_end);
    //     real_t delta_u = (u_end - u)/ N;

    //     for (int i = 0; i < N + 1; i++) {
    //         fc = exp(u + i*delta_u);
    //         real_t PLoss = ComputePowerLoss(fc, fc_mu);
    //         data_file << fc << ";" << PLoss << std::endl; 
            
            
    //     }
    // }


    // real_t fc = 500e6;
    // real_t fc_mu = 1.8e6;
    // ComputePowerLoss(fc, fc_mu);











    std::string name = "./data";
    name += std::to_string(0) + ".csv";

    // ofstream to those files
    std::ofstream data_file(name);
    data_file << "fc;Ploss;flux\n";
    real_t fc_0 = 100e3;
    real_t fc = fc_0;
    real_t fc_end = 2e6;
    int N = 24; // N+1 points

    real_t fc_mu = 1.8e6;

    // Change of variable to obtain points linearly spaced on a logarithmic scale
    real_t u = log(fc_0);
    real_t u_end = log(fc_end);
    real_t delta_u = (u_end - u)/ N;

    for (int i = 0; i < N + 1; i++) {
        fc = exp(u + i*delta_u);                            // Frequency for the simulation
        // real_t PLoss, flux;                                 // parameters computed
        real_t imax = 0;
        real_t Ploss = ComputePowerLoss(fc, fc_mu);         // computation of the parameters
        data_file << fc << ";" << Ploss << std::endl; 
    }
}