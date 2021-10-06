#include "types.h"

#include <stdio.h>

void mass_damper_spring(Vector * p_xdot, const Vector * p_x){

}

void runge_kutta_45(Plant plant, Vector * p_x, Vector * p_xtemp){
    static const double butcher_a[][6] = {
        {1.0/5},
        {3.0/40, 9.0/40},
        {44.0/45, -56.0/15, 32.0/9},
        {19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729},
        {9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656},
        {35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84}
    };

    static const double butcher_b[] = {
        35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84
    };

    static const double butcher_c[] = {
        1.0/5, 3.0/10, 4.0/5, 8.0/9, 1.0, 1.0
    };

    static double Vector


    for(int i = 0; i < p_x->size; i++){
        p_xtemp->data[i] = p_x->data[i];
    }
}

int main(){
    return 0;
}
