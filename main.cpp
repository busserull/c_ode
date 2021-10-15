#include <iostream>
#include <memory>

// ------ Vec

class Vec {
public:
    Vec(int dim);
    double & operator[](int dim);
private:
    int m_dim;
    std::unique_ptr<double[]> m_data;
};

Vec::Vec(int dim){
    m_dim = dim;
    m_data = std::unique_ptr<double[]>(new double[dim]);

    for(int i = 0; i < dim; i++){
        m_data[i] = 0.0;
    }
}

double & Vec::operator[](int dim){
    return m_data[dim];
}



int main(){
    Vec a(2);

    std::cout << a[0] << " " << a[1] << std::endl;

    a[0] = 12.0;
    a[1] = 6.32;

    std::cout << a[0] << " " << a[1] << std::endl;

    return 0;
}
