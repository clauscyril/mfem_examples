#ifndef FERRITE
#define FERRITE
#include <mfem.hpp>

#ifndef M_PI 
#define M_PI 3.14159265358979323846 
#endif

using namespace mfem;
struct Ferrite
{
    std::string name;
    real_t rho;
    real_t sigma;
    real_t eps;
    real_t mu; 
    Ferrite(const char* name_, real_t rho_, real_t sigma_, real_t eps_, real_t mu_r_) // IMPORTANT IT IS THE RELATIVE PERMITTIVITY
    : name(name_),  rho(rho_), sigma(sigma_), eps(eps_), mu(4e-7*M_PI*mu_r_){};       // THAT IS PASSED AS AN ARGUMENT 
};
#endif