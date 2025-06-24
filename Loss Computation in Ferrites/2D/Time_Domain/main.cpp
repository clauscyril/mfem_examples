#include <iostream>
#include <complex>
#include "mfem.hpp"
#include "solver_TD.hpp"


using namespace mfem;


int main() {
    const char *path = "../../../mesh/square.msh";

    Mesh *mesh = new Mesh(path, 1, 1); 
    mesh->UniformRefinement();
    // mesh->UniformRefinement();
    // mesh->UniformRefinement();

    Ferrite N30("N30", 5.98e-2, 4.44e-1, 2.48e-6, 4300);
    Ferrite N87("N87", 4.24e-2, 1.48e-1, 2.68e-6, 2200);
    Ferrite T38("T38", 4.04e-2, 1.07e1, 8.06e-6, 10000);

    real_t height = 7.59e-3;
    real_t w = 5.3e-3;
    real_t b_peak = 10e-3; // Average value of b
    real_t phi_peak = b_peak*w*height;
    
    // real_t I_rms = 0.082348/sqrt(2);
    real_t f = 500e3;
    real_t I_rms = 0.082/sqrt(2);
    int nb_period = 5;
    int nb_iter = 1000 * nb_period;
    real_t Ts = nb_period/f/nb_iter;

    real_t t_f = nb_period/f;

    Ferrite ferrite = N30;
    real_t tau = 1/(2*M_PI*1.8e6);

    real_t rho = ferrite.rho;
    real_t sigma = ferrite.sigma;
    real_t eps = ferrite.eps;
    real_t mu = ferrite.mu; 

    real_t A1 = (rho * (sigma * Ts + 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A2 = (rho * (sigma * Ts - 2*eps) + Ts)/(2*eps + Ts*sigma);
    real_t A3 = -(Ts*sigma -2*eps) / (2*eps + Ts*sigma);

    real_t B1 = 2*mu / (Ts + 2*tau);
    real_t B2 = -(2*mu) / (Ts + 2*tau);
    real_t B3 = -(Ts - 2*tau) / (Ts + 2*tau);

    real_t C1 = Ts*mu/(Ts+2*tau);
    real_t C2 = Ts*mu/(Ts+2*tau);
    real_t C3 = -(Ts-2*tau)/(Ts+2*tau);

    // std::cout << "A1 = " << A1 << std::endl;
    // std::cout << "A2 = " << A2 << std::endl;
    // std::cout << "A3 = " << A3 << std::endl;
    // std::cout << "B1 = " << B1 << std::endl;
    // std::cout << "B2 = " << B2 << std::endl;
    // std::cout << "B3 = " << B3 << std::endl;
    // std::cout << "C1 = " << C1 << std::endl;
    // std::cout << "C2 = " << C2 << std::endl;
    // std::cout << "C3 = " << C3 << std::endl;
    // 
    auto NI_sine_func = [&](real_t t) 
    {
        return I_rms * std::sqrt(2) * std::sin(2 * M_PI * f * t);
    };

    auto NI_square_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return I_rms * std::sqrt(2)/2;
        else
            return -I_rms * std::sqrt(2)/2.;
    };

    
    auto NI_square_pos_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return I_rms * std::sqrt(2);
        else
            return 0.;
    };


    auto NI_saw_pos_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        return I_rms * std::sqrt(2) * rest * Ts;
    };  

    auto NI_saw_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return I_rms * std::sqrt(2) * rest * Ts;
        else
            return I_rms * std::sqrt(2) * (nb_iter/nb_period - rest) * Ts;
    };
    auto NI_zero = [&](real_t t) 
    {
        return 0;
    };


    auto phi_sine_func = [&](real_t t) 
    {
        return phi_peak * std::sin(2 * M_PI * f * t);
    };

    auto phiH_sine_func = [&](real_t t) 
    {
        return 7.7805e-5 * std::sin(2 * M_PI * f * t);
    };



    std::complex<real_t> j(0,1);
    real_t omega = 2*M_PI*f;           // Working frequency
    const std::complex<real_t> mu_eq = mu / ((real_t)1. + tau * j * omega);
    real_t phase = std::abs(mu_eq);
    auto phi_sine_func_phased = [&](real_t t) 
    {
        return phi_peak * std::sin(2 * M_PI * f * t - phase);
    };

    auto phi_square_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return phi_peak /2;
        else
            return -phi_peak/2.;
    };

    auto phi_saw_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return phi_peak * rest * Ts;
        else
            return phi_peak* (nb_iter/nb_period - rest) * Ts;
    };

    auto phiH_saw_func = [&](real_t t) 
    {
        int iter = int(t/Ts);
        int rest = iter%((int)(nb_iter/nb_period));

        if (rest < nb_iter/nb_period/2)
            return 7.7805e-5 * rest * Ts;
        else
            return 7.7805e-5* (nb_iter/nb_period - rest) * Ts;
    };

    // TD_sim(mesh, NI_sine_func, t_f, nb_iter, N30, false);
    TD_sim_by_flux(mesh, phi_saw_func, t_f, nb_iter, N30, false);

    // *********** TEST ***********
    // real_t phi_n = 0;
    // real_t phi_nm1 = 0;
    // real_t phiH_n = 0;
    // real_t phiH_nm1 = 0;

    // real_t t = 0;
    // std::string name_test = "./test.csv";   
    // std::ofstream data_file_test(name_test);
    // data_file_test << "t;phi;phiH" << std::endl;
    // for (int step = 0; step < nb_iter; step++)
    // {   
    //     t += Ts;
    //     phi_n = phi_sine_func(t);
    //     phiH_n = 1/C1*(phi_n - C2*phiH_nm1 - C3*phi_nm1);
    //     phiH_nm1 = phiH_n;
    //     phi_nm1 = phi_n;
    //     data_file_test << t << ";" << phi_n  << ";" << phiH_n << std::endl; 
    // }
    return 0;
}


