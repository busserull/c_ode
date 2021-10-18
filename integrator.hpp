#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP
#include <memory>
#include "vec.hpp"

enum class IntegratorMethod {
    explicit_euler,
    midpoint,
    rk4,
    dormand_prince,
};

class Integrator {
public:
    Integrator(IntegratorMethod method, int plant_dimension);
    void step(Vec &x_out, double t, const Vec &x, const Vec &u, double h);
private:
    int m_steps;
    std::unique_ptr<double[]> m_butcher_a;
    std::unique_ptr<double[]> m_butcher_b;
    std::unique_ptr<double[]> m_butcher_c;
    std::unique_ptr<double[]> m_butcher_e;
    std::unique_ptr<Vec[]> m_k;
    Vec m_stage;
    Vec m_temp;
};

#endif
