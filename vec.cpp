#include "vec.hpp"

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
