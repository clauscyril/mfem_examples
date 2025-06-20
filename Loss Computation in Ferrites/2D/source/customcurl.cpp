#include "../headers/customcurl.hpp"

using namespace mfem;
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

// if(0)
//     {
//         V[0] = -grad_H[1];                  // - dH/dz
//         V[1] = grad_H[0] + H_val /(x[0]);   // dH/dr + H/r
//     }
//     {
//         int order = 1;
//         int dim = 2;
//         FiniteElementSpace fespace_h1 = const_cast<FiniteElementSpace> (GridFunc->FESpace());
//         ND_FECollection fe_coll_nd(order, dim);
//         // FiniteElementSpace fespace_nd(GridFunc->FESpace()->GetMesh(), &fe_coll_nd);
//         FiniteElementSpace fespace_nd(fespace_h1->GetMesh(), &fe_coll_nd);
//         GridFunction dx(&fespace_nd);
//         DiscreteLinearOperator *grad = new DiscreteLinearOperator(fespace_h1, &fespace_nd);
//         grad->AddDomainInterpolator(new GradientInterpolator);
//         grad->Assemble();
//         grad->Finalize();
//         grad->Mult( *(GridFunc), dx);    
//         V = dx.GetTrueVector();        
//     }