#include <iostream>
#include <fstream>
// #include <algorithm>

#include <mfem.hpp>

using namespace mfem;
using namespace std;


real_t q_func(const Vector &x)
{
    if (abs(x(0)) < 0.2 and abs(x(1)) < 0.2) 
        return 1.;
    return 0.;
}

class DG_Solver : public Solver
{
private:
   SparseMatrix &M, &K, A;
   GMRESSolver linear_solver;
   BlockILU prec;
   real_t dt;
public:
   DG_Solver(SparseMatrix &M_, SparseMatrix &K_, const FiniteElementSpace &fes)
      : M(M_),
        K(K_),
        prec(fes.GetTypicalFE()->GetDof(),
             BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
        dt(-1.0)
   {
      linear_solver.iterative_mode = false;
      linear_solver.SetRelTol(1e-9);
      linear_solver.SetAbsTol(0.0);
      linear_solver.SetMaxIter(100);
      linear_solver.SetPrintLevel(0);
      linear_solver.SetPreconditioner(prec);
   }

   void SetTimeStep(real_t dt_)
    {
        if (dt_ != dt)
        {
            dt = dt_;
            // Form operator A = M - dt*K
            A = K;
            A *= -dt;
            A += M;

            // this will also call SetOperator on the preconditioner
            linear_solver.SetOperator(A);
        }
    }

   void SetOperator(const Operator &op) override
    {
        linear_solver.SetOperator(op);
    }

   void Mult(const Vector &x, Vector &y) const override
   {
        linear_solver.Mult(x, y);
   }
};

/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
    class FE_Evolution : public TimeDependentOperator
    {
    private:
       BilinearForm &M, &K;
       const Vector &b;
       Solver *M_prec;
       CGSolver M_solver;
       DG_Solver *dg_solver;
    
       mutable Vector z;
    
    public:
       FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_);
    
       void Mult(const Vector &x, Vector &y) const override;
       void ImplicitSolve(const real_t dt, const Vector &x, Vector &k) override;
    
       ~FE_Evolution() override;
    };
    

int main() {

    int vis_steps = 5;
    bool visualization = true;
    int precision = 8;

    real_t t_final = 5.;
    real_t dt = 0.05;

    real_t rho_ = 1.;
    real_t c_ = 1.;

    real_t k_ = 1.;

    unique_ptr<ODESolver> ode_solver = ODESolver::Select(4); // Runge kuta

    int n;       
    int order = 1;  

    const char *path = "../Disque2.msh";
    Mesh mesh(path, 1, 1);

    int ne = mesh.GetNE();
    int dim = mesh.Dimension();
    int spaceDim = mesh.SpaceDimension();

    cout << "Nombre d'elements : " << ne << endl;
    cout << "Dim : " << dim << "\nSpaceDim : " << spaceDim << endl;


    FiniteElementCollection *fec = new H1_FECollection(order, dim);    
    FiniteElementSpace *fespace = new FiniteElementSpace(&mesh, fec);


    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Liste des dofs de bord
    cout << "Taille bdr_attribute : " << mesh.bdr_attributes.Size() << endl;


    GridFunction T(fespace);


    FunctionCoefficient q_function(q_func);
    LinearForm q = LinearForm(fespace);
    q.AddDomainIntegrator(new DomainLFIntegrator(q_function));


    BilinearForm m = BilinearForm(fespace);
    BilinearForm k = BilinearForm(fespace);

    ConstantCoefficient alpha(rho_*c_);
    m.AddDomainIntegrator(new MassIntegrator(alpha));

    ConstantCoefficient k_coef(k_);
    k.AddDomainIntegrator(new DiffusionIntegrator(k_coef));

    k.Assemble();
    m.Assemble();

    k.Finalize(0);
    m.Finalize(0);

    q.Assemble(); 


    socketstream sout;

    char vishost[] = "localhost";
    int  visport   = 19916;
    sout.open(vishost, visport);
    if (!sout)
    {
        cout << "Unable to connect to GLVis server at "
            << vishost << ':' << visport << endl;
        visualization = false;
        cout << "GLVis visualization disabled.\n";
    }
    else
    {
        sout.precision(precision);
        sout << "solution\n" << mesh << T;
        sout << "pause\n";
        sout << flush;
        cout << "GLVis visualization paused."
            << " Press space (in the GLVis window) to resume it.\n";
    }
    



    // Solver
    FE_Evolution adv(m, k, q);

    real_t t = 0.0;
    adv.SetTime(t);
    ode_solver->Init(adv);



    bool done = false;
    for (int ti = 0; !done; )
    {
        real_t dt_real = min(dt, t_final - t);
        ode_solver->Step(T, t, dt_real);
        ti++;

        done = (t >= t_final - 1e-8*dt);

        if (done || ti % vis_steps == 0)
        {
            cout << "time step: " << ti << ", time: " << t << endl;
            
            if (visualization)
            {
                sout << "solution\n" << mesh << T << flush;
            }
        
        }
    }

    return 0;
}





// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(BilinearForm &M_, BilinearForm &K_, const Vector &b_)
   : TimeDependentOperator(M_.FESpace()->GetTrueVSize()),
     M(M_), K(K_), b(b_), z(height)
{
   Array<int> ess_tdof_list;
   if (M.GetAssemblyLevel() == AssemblyLevel::LEGACY)
   {
      M_prec = new DSmoother(M.SpMat());
      M_solver.SetOperator(M.SpMat());
      dg_solver = new DG_Solver(M.SpMat(), K.SpMat(), *M.FESpace());
   }
   else
   {
      M_prec = new OperatorJacobiSmoother(M, ess_tdof_list);
      M_solver.SetOperator(M);
      dg_solver = NULL;
   }
   M_solver.SetPreconditioner(*M_prec);
   M_solver.iterative_mode = false;
   M_solver.SetRelTol(1e-9);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
   // y = M^{-1} (K x + b)
   K.Mult(x, z);
   z += b;
   M_solver.Mult(z, y);
}

void FE_Evolution::ImplicitSolve(const real_t dt, const Vector &x, Vector &k)
{
   MFEM_VERIFY(dg_solver != NULL,
               "Implicit time integration is not supported with partial assembly");
   K.Mult(x, z);
   z += b;
   dg_solver->SetTimeStep(dt);
   dg_solver->Mult(z, k);
}

FE_Evolution::~FE_Evolution()
{
   delete M_prec;
   delete dg_solver;
}