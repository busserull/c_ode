#include <iostream>
#include <memory>

// ------ Vec

class Vec {
public:
    Vec();
    Vec(int dim);

    int dim() const;
    void null();

    void operator=(const Vec &rhs);
    void operator+=(const Vec &rhs);
    void operator*=(double rhs);

    double & operator[](int dim);
    double operator[](int dim) const;
private:
    int m_dim;
    std::unique_ptr<double[]> m_data;
};

Vec::Vec(){
    m_dim = 0;
    m_data = nullptr;
}

Vec::Vec(int dim){
    m_dim = dim;
    m_data = std::unique_ptr<double[]>(new double[dim]);

    for(int i = 0; i < dim; i++){
        m_data[i] = 0.0;
    }
}

int Vec::dim() const {
    return m_dim;
}

void Vec::null(){
    for(int i = 0; i < m_dim; i++){
        m_data[i] = 0.0;
    }
}

void Vec::operator=(const Vec &rhs){
    if(m_dim != rhs.m_dim){
        m_dim = rhs.m_dim;
        m_data = std::unique_ptr<double[]>(new double[m_dim]);
    }

    // Assert this->dim == rhs.dim
    for(int i = 0; i < m_dim; i++){
        m_data[i] = rhs.m_data[i];
    }
}

void Vec::operator+=(const Vec &rhs){
    // Assert this->dim == rhs.dim
    for(int i = 0; i < m_dim; i++){
        m_data[i] += rhs.m_data[i];
    }
}

void Vec::operator*=(double rhs){
    for(int i = 0; i < m_dim; i++){
        m_data[i] *= rhs;
    }
}

double & Vec::operator[](int dim){
    return m_data[dim];
}

double Vec::operator[](int dim) const {
    return m_data[dim];
}


// ------ Integrator

enum class IntegratorMethod {
    explicit_euler,
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

Integrator::Integrator(IntegratorMethod method, int plant_dimension){
    switch(method){
        case IntegratorMethod::explicit_euler:
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
    x_dot[0] = x[1];
    x_dot[1] = -1.0 * x[0] - 0.5 * x[1];
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



int main(){
    Vec x(2);
    x[0] = 1.0;
    x[1] = 0.0;

    Integrator s(IntegratorMethod::rk4, 2);

    std::cout << "x = [..." << std::endl;
    for(int i = 0; i < 300; i++){
        std::cout << x[0] << ",..." << std::endl;
        s.step(x, 0.0, x, x, 0.1);
    }
    std::cout << "];" << std::endl;

    std::cout << "plot(x);" << std::endl;

    return 0;
}
