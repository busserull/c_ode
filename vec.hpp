#ifndef VEC_HPP
#define VEC_HPP
#include <memory>

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

#endif
