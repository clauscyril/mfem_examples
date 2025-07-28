#ifndef CURLCOEFF
#define CURLCOEFF
#include "mfem.hpp"
using namespace mfem;


// custom curl of an axisymetric problem, with H = H_theta . e_theta and curl(H) = -dH_theta/dz . e_r + (dH_theta/dr + H_tehta/r) e_z
// This class is inspired from GradientGridFunctionCoefficient
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
#endif