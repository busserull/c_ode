#include <iostream>
#include "vec.hpp"
#include "integrator.hpp"




void integrate_plot(Integrator &s, double end_time, double step_size){
    double time = 0.0;

    Vec x(2);
    x[0] = 1.0;
    x[1] = 0.0;

    std::cout << "x = [..." << std::endl;
    while(time < end_time){
        std::cout << x[0] << ",..." << std::endl;
        s.step(x, 0.0, x, x, step_size);
        time += step_size;
    }
    std::cout << "]; plot(x);" << std::endl;
}


int main(){
    double end_time = 5.0;
    double step_size = 1.0;

    std::cout << "figure(); hold on;" << std::endl;

    Integrator s_explicit(IntegratorMethod::explicit_euler, 1);
    Integrator s_midpoint(IntegratorMethod::midpoint, 1);
    Integrator s_rk4(IntegratorMethod::rk4, 1);
    Integrator s_dormand_price(IntegratorMethod::dormand_prince, 1);

    integrate_plot(s_explicit, end_time, step_size);
    integrate_plot(s_midpoint, end_time, step_size);
    integrate_plot(s_rk4, end_time, step_size);
    integrate_plot(s_dormand_price, end_time, step_size);

    std::cout << "legend('Explicit', 'Midpoint', 'RK4', 'Dormand Prince');" << std::endl;

    return 0;
}
