#include "integrator.hpp"

Integrator::Integrator(IntegratorMethod method, int plant_dimension){
    switch(method){
        case IntegratorMethod::explicit_euler:
            m_steps = 1;
            m_butcher_a = nullptr;
            m_butcher_b = std::unique_ptr<double[]>(new double[1]);
            m_butcher_c = std::unique_ptr<double[]>(new double[1]);
            m_butcher_e = nullptr;

            m_butcher_b[0] = 1.0;

            m_butcher_c[0] = 0.0;
            break;

        case IntegratorMethod::midpoint:
            m_steps = 2;
            m_butcher_a = std::unique_ptr<double[]>(new double[1]);
            m_butcher_b = std::unique_ptr<double[]>(new double[2]);
            m_butcher_c = std::unique_ptr<double[]>(new double[2]);
            m_butcher_e = nullptr;

            m_butcher_a[0] = 1.0/2.0;

            m_butcher_b[0] = 0.0;
            m_butcher_b[1] = 1.0;

            m_butcher_c[0] = 0.0;
            m_butcher_c[1] = 1.0/2.0;
            break;

        case IntegratorMethod::rk4:
            m_steps = 4;
            m_butcher_a = std::unique_ptr<double[]>(new double[6]);
            m_butcher_b = std::unique_ptr<double[]>(new double[4]);
            m_butcher_c = std::unique_ptr<double[]>(new double[4]);
            m_butcher_e = nullptr;

            m_butcher_a[0] = 0.5;
            m_butcher_a[1] = 0.0;
            m_butcher_a[2] = 0.5;
            m_butcher_a[3] = 0.0;
            m_butcher_a[4] = 0.0;
            m_butcher_a[5] = 1.0;

            m_butcher_b[0] = 1.0/6.0;
            m_butcher_b[1] = 1.0/3.0;
            m_butcher_b[2] = 1.0/3.0;
            m_butcher_b[3] = 1.0/6.0;

            m_butcher_c[0] = 0.0;
            m_butcher_c[1] = 0.5;
            m_butcher_c[2] = 0.5;
            m_butcher_c[3] = 1.0;
            break;

        case IntegratorMethod::dormand_prince:
            m_steps = 7;
            m_butcher_a = std::unique_ptr<double[]>(new double[21]);
            m_butcher_b = std::unique_ptr<double[]>(new double[7]);
            m_butcher_c = std::unique_ptr<double[]>(new double[7]);
            m_butcher_e = std::unique_ptr<double[]>(new double[7]);

            m_butcher_a[0] = 1.0/5.0;
            m_butcher_a[1] = 3.0/40.0;
            m_butcher_a[2] = 9.0/40.0;
            m_butcher_a[3] = 44.0/45.0;
            m_butcher_a[4] = -56.0/15.0;
            m_butcher_a[5] = 32.0/9.0;
            m_butcher_a[6] = 19372.0/6561.0;
            m_butcher_a[7] = -25360.0/2187.0;
            m_butcher_a[8] = 64448.0/6561.0;
            m_butcher_a[9] = -212.0/729.0;
            m_butcher_a[10] = 9017.0/3168.0;
            m_butcher_a[11] = -355.0/33.0;
            m_butcher_a[12] = 46732.0/5247.0;
            m_butcher_a[13] = 49.0/176.0;
            m_butcher_a[14] = -5103.0/18656.0;
            m_butcher_a[15] = 35.0/384.0;
            m_butcher_a[16] = 0.0;
            m_butcher_a[17] = 500.0/1113.0;
            m_butcher_a[18] = 125.0/192.0;
            m_butcher_a[19] = -2187.0/6784.0;
            m_butcher_a[20] = 11.0/84.0;

            m_butcher_b[0] = 35.0/384.0;
            m_butcher_b[1] = 0.0;
            m_butcher_b[2] = 500.0/1113.0;
            m_butcher_b[3] = 125.0/192.0;
            m_butcher_b[4] = -2187.0/6784.0;
            m_butcher_b[5] = 11.0/84.0;
            m_butcher_b[6] = 0.0;

            m_butcher_c[0] = 0.0;
            m_butcher_c[1] = 1.0/5.0;
            m_butcher_c[2] = 3.0/10.0;
            m_butcher_c[3] = 4.0/5.0;
            m_butcher_c[4] = 8.0/9.0;
            m_butcher_c[5] = 1.0;
            m_butcher_c[6] = 1.0;

            m_butcher_e[0] = 5179.0/57600.0;
            m_butcher_e[1] = 0.0;
            m_butcher_e[2] = 7571.0/16695.0;
            m_butcher_e[3] = 393.0/640.0;
            m_butcher_e[4] = -92097.0/339200.0;
            m_butcher_e[5] = 187.0/2100.0;
            m_butcher_e[6] = 1.0/40.0;
            break;
    }

    m_k = std::unique_ptr<Vec[]>(new Vec[m_steps]);

    for(int i = 0; i < m_steps; i++){
        m_k[i] = Vec(plant_dimension);
    }

    m_stage = Vec(plant_dimension);
    m_temp = Vec(plant_dimension);
}

void plant(Vec &x_dot, double _t, const Vec &x, const Vec &u){
    x_dot[0] = x[0];
    /* x_dot[0] = x[1]; */
    /* x_dot[1] = -1.0 * x[0] - 0.5 * x[1]; */
}

void Integrator::step(Vec &x_out, double t, const Vec &x, const Vec &u, double h){
    int butcher_a_index = 0;

    for(int step = 0; step < m_steps; step++){
        m_stage.null();

        for(int i = 0; i < step; i++){
            m_temp = m_k[i];
            m_temp *= m_butcher_a[butcher_a_index];
            m_stage += m_temp;

            butcher_a_index++;
        }

        m_stage *= h;
        m_stage += x;

        // Perform plant equation evaluation
        // Put into k for step
        plant(m_k[step], t, x, u);
    }

    m_stage.null();

    for(int i = 0; i < m_steps; i++){
        m_temp = m_k[i];
        m_temp *= m_butcher_b[i];
        m_stage += m_temp;
    }

    m_stage *= h;
    m_stage += x;

    x_out = m_stage;
}
