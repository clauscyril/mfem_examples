#include <iostream>
#include <complex>
#include "mfem.hpp"


#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

using namespace mfem;

// NI parameters
real_t I_max = 0.08;
real_t f = 100e3;


// geometric parameters
real_t Ri = 9.6e-3/2.0;
real_t height = 7.59e-3;
real_t w = 5.3e-3;
real_t Rout = Ri + w;
real_t Rm = (Rout - Ri) / 2;

// ****** Ferrite parameters ******

// Ferrite N30
real_t rho = 5.98e-2;
real_t sigma = 4.44e-1;
real_t eps = 2.48e-6;
real_t mu = 4300.0 * 4e-7 * M_PI;

real_t tau = 1/(2 * M_PI * 1.8e6);

int num_steps = 10000;
real_t t = 0.0;
real_t Ts = 1./f/1000;



real_t A3 = sigma * Ts + 2*eps;
real_t A4 = sigma * Ts - 2*eps;
real_t A1 = rho * A3 + Ts;
real_t A2 = rho * A4 + Ts;

// real_t B1 = 1;
// real_t B2 = 1;
// real_t B3 = 1;
// real_t B4 = 1;

real_t C1 = 2*mu;
real_t C2 = Ts + 2*tau;
real_t C3 = Ts - 2*tau;



real_t NI = 0;

real_t bdr_func(const Vector &x, real_t &t);

real_t r_coeff_func(const Vector &x);

real_t inv_r_square_func(const Vector &x);


class CurlCustomCoefficient : public VectorCoefficient
{
protected:
   const GridFunction *GridFunc;

public:
   CurlCustomCoefficient(const GridFunction *gf);
   void SetGridFunction(const GridFunction *gf);
   const GridFunction * GetGridFunction() const { return GridFunc; }

   void Eval(Vector &V, ElementTransformation &T,
             const IntegrationPoint &ip) override;

   virtual ~CurlCustomCoefficient() { }
};




int main(){
    const char *path = "../../../mesh/square.msh";             // Path to the mesh

    Mesh *mesh = new Mesh(path, 1, 1);
    // mesh->UniformRefinement();

    int order = 2;                      // elements order : 2
    int dim = mesh->Dimension();        // Mesh dimension : 2

    FiniteElementCollection *fec = new H1_FECollection(order, dim); 
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

    FiniteElementSpace *fespace_E = new FiniteElementSpace(mesh, fec, dim);

    GridFunction Hn(fespace);     // H^n
    GridFunction Hnm1(fespace);   // H^{n-1}
    GridFunction db_dt(fespace);
    GridFunction db_dtm1(fespace);

    GridFunction En(fespace_E);
    GridFunction Enm1(fespace_E);

    GridFunction Jn(fespace_E);
    GridFunction Jnm1(fespace_E);
    // CurlCustomCoefficient *curl_H_coeff = new CurlCustomCoefficient(&Hn);

    En = 0;
    Jn = 0;
    Enm1 = 0;
    Enm1 = 0;

    Hnm1 = 0;  // Initial values
    db_dt = 0;

    // Definition of the Boundary conditions  :
    Array<int> ess_tdof_list;
    Array<int> dir_bdr(mesh->bdr_attributes.Max());
    dir_bdr = 1; // All the borders have boundary conditions
    fespace->GetEssentialTrueDofs(dir_bdr, ess_tdof_list);

    FunctionCoefficient bdr_coeff([&](const Vector &x) { return bdr_func(x, t); });
    Hnm1.ProjectBdrCoefficient(bdr_coeff, dir_bdr);

    // **** Coefficients for bilinear forms
    FunctionCoefficient r_coeff(r_coeff_func);
    FunctionCoefficient inv_r_square_coeff(inv_r_square_func);

    ConstantCoefficient A1C2_coeff(A1*C2);
    ProductCoefficient A1C2_r(A1C2_coeff, r_coeff);

    ConstantCoefficient A3C1(C1*A3);
    ProductCoefficient A3C1_r(A3C1, r_coeff);

    ProductCoefficient A1C2_r_inv(A1C2_r, inv_r_square_coeff);


    ConstantCoefficient A2C2_coeff(-A2*C2);
    ProductCoefficient A2C2_r(A2C2_coeff, r_coeff);

    ConstantCoefficient A3C1m(-C1*A3);
    ProductCoefficient A3C1m_r(A3C1m, r_coeff);

    ProductCoefficient A2C2_r_inv(A2C2_r, inv_r_square_coeff);

    ConstantCoefficient k3(C3*A3 - A4*C2);
    ProductCoefficient k3_r(k3, r_coeff);


    BilinearForm r1(fespace); 
    r1.AddDomainIntegrator(new DiffusionIntegrator(A1C2_r));
    r1.AddDomainIntegrator(new MassIntegrator(A3C1_r));
    r1.AddDomainIntegrator(new MassIntegrator(A1C2_r_inv));
    r1.Assemble();

    BilinearForm r2(fespace); 
    r2.AddDomainIntegrator(new DiffusionIntegrator(A2C2_r));
    r2.AddDomainIntegrator(new MassIntegrator(A3C1m_r));
    r2.AddDomainIntegrator(new MassIntegrator(A2C2_r_inv));
    r2.Assemble();


    // Initialisation à 0
    // GridFunctionCoefficient db_dt_coeff(&db_dt); 
    // LinearForm *m = new LinearForm(fespace);
    // m->AddDomainIntegrator(new DomainLFIntegrator(k3_r));
    // m->Assemble();
    BilinearForm m(fespace);
    m.AddDomainIntegrator(new MassIntegrator(k3_r));
    m.Assemble();


    r1.Finalize();
    r2.Finalize();
    m.Finalize();



    SparseMatrix &A = r1.SpMat();
    Vector rhs(fespace->GetTrueVSize());
    Vector R2Hnm1(fespace->GetTrueVSize());
    Vector MdBdt(fespace->GetTrueVSize());

    Vector x(fespace->GetTrueVSize()); // Pour stocker solution complète

    CGSolver solver;
    solver.SetOperator(A);
    solver.SetRelTol(1e-12);
    solver.SetMaxIter(5000);
    solver.SetPrintLevel(0);


    socketstream sout;
    socketstream sout_e;
    socketstream sout_j;
    char vishost[] = "localhost";
    int visport = 19916;
    sout.open(vishost, visport);
    sout_e.open(vishost, visport);
    sout_j.open(vishost, visport);

    sout.precision(8);
    sout << "solution\n" << *mesh << Hn << "\nkeys j\n" << std::flush;
    sout_e.precision(8);
    sout_e << "solution\n" << *mesh << En << "\nkeys j\n" << std::flush;
    sout_j.precision(8);
    sout_j << "solution\n" << *mesh << Jn << "\nkeys j\n" << std::flush;


    int vis_steps = 10;
    for (int step = 0; step < num_steps; step++)
    {
        t += Ts;
        // std::cout << "H(0) = " << Hn(0) <<  std::endl;
        // std::cout << "t = " << t << std::endl;
        std::cout << "NI/(2*pi*Ri) = " << I_max * sqrt(2)*sin(2*M_PI*f * t) /(2*M_PI*Ri) << std::endl;
        std::cout << "A3 = " << A3 << std::endl;
        // Mettre à jour BDConditions avec ProjectCoefficient
        // ....
        // calcul bdr_coeff
       // Mettre à jour le coefficient de bord avec le temps courant

        bdr_coeff.SetTime(t);

        // Appliquer la condition de Dirichlet à Hn sur les bords
        Hn.ProjectBdrCoefficient(bdr_coeff, dir_bdr);

        // Mult M * dBdt → MdBdt
        m.Mult(db_dt, MdBdt);

        // Mult R2 * Hnm1 → R2Hnm1
        r2.Mult(Hnm1, R2Hnm1);

        // RHS = R2Hnm1 + MdBdt
        rhs = R2Hnm1;
        rhs += MdBdt;
        
        SparseMatrix A_sys;
        Vector X, B;

        r1.FormLinearSystem(ess_tdof_list, Hn, rhs, A_sys, X, B);

        // Résolution du système linéaire : A_sys * X = B
        solver.SetOperator(A_sys);
        solver.Mult(B, X);

        // Reconstruire Hn à partir de X
        r1.RecoverFEMSolution(X, rhs, Hn);

        // Mise à jour des champs dérivés (db_dt, Hnm1)
        db_dtm1 *= C3 / C2;

        Hn *= C1 / C3;
        Hnm1 *= C1 / C2;

        db_dt = 0;
        db_dt += Hn;
        db_dt -= Hnm1;
        db_dt -= db_dtm1;

        Hnm1 = Hn;
        Hnm1 *= C3 / C1;
        Hn *= C3/C1;




        CurlCustomCoefficient J_coeff(&Hn);
        Jn.ProjectCoefficient(J_coeff);

        Enm1 *= A4;
        Jnm1 *= A2;
        Jn *= A1;

        En += Jn;
        En += Jnm1;
        En -= Enm1;

        En *= 1./A3;
        Jn *= 1./A1;

        Enm1 = En;
        Jnm1 = Jn;

        bool visualization = true;

        

        // Visualisation avec GLVis
        if (step % vis_steps == 0)
        {
            std::cout << "Time step: " << step << ", Time: " << t << std::endl;

            if (visualization)
            {
                sout.precision(8);
                sout << "solution\n"
                     << *mesh << Hn
                     << "window_title 'Champ H'"
                     << std::flush;
                sout_e << "solution\n" << *mesh << En <<  "window_title 'Champ E'" << std::flush;
                sout_j << "solution\n" << *mesh << Jn <<  "window_title 'Courant J'" << std::flush;
                // std::cin.get(); // Décommenter pour pause manuelle
            }
        }
        
    }

}


// Boundary conditions function
real_t bdr_func(const Vector &x, real_t &t)
{
    real_t r = x[0];
    NI = I_max * sqrt(2)*sin(2*M_PI*f * t);
    return NI /(2 * M_PI * r); 
}

// Function used for the integration because of the cylindrical coordinates
real_t r_coeff_func(const Vector &x){
    return x[0];  
}

// Function for the correcting factor rho/r² 
real_t inv_r_square_func(const Vector &x){
    return (real_t)1./pow(x[0],2);
}


// Constructor of custom curl of an axisymetric problem 
CurlCustomCoefficient::CurlCustomCoefficient (const GridFunction *gf)  // takes only a gridFunction as argument
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
    Vector x;           // Vector of the coordinates (r, z)
    T.Transform(ip, x); // Get the coordinates of the point ip in vector x, because ip has its coordinates in the element's referential                

    real_t H_val = GridFunc->GetValue(T, ip);   // Get value of H in indicated coordinates

    Vector grad_H(2);     // Gradient vector            
    GridFunc->GetGradient(T, grad_H); // Evaluating the gradient value at the integration point ip
                                      // In this case, ip is already linked to the ElementTransformation

    V[0] = -grad_H[1];                  // - dH/dz
    V[1] = grad_H[0] + H_val /(x[0]);   // dH/dr + H/r

} 