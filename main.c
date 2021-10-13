#include "types.h"
#include "integrator.h"

#include <stdio.h>
#include <stdlib.h>

void integrator_step(Integrator *, Vector *, double, double, const Vector *, const Vector *);

void mass_damper_spring(
    void * p_params, 
    Vector * p_xdot, 
    double t, 
    const Vector * p_x, 
    const Vector * p_u
){
    p_xdot->data[0] = p_x->data[1];
    p_xdot->data[1] = -1 * p_x->data[0] - 0.1 * p_x->data[1];
}

/* void mass_damper_spring(Vector * p_xdot, const Vector * p_x){ */

/* } */

int main(){
    Plant mds;
    mds.dim = 2;
    mds.p_params = NULL;
    mds.xdot = mass_damper_spring;

    Integrator s = integrator_new(&mds, INTEGRATOR_METHOD_RK4);

    Vector x;
    x.dim = 2;
    x.data = (double *)malloc(2 * sizeof(double));
    x.data[0] = 1.0;
    x.data[1] = 0.0;

    double t = 0;
    double h = 0.1;

    while(t < 2.0){
        printf("%f\n", x.data[0]);
        integrator_step(&s, &x, h, t, &x, NULL);
        t += h;
    }

    free(x.data);

    integrator_delete(&s);

    return 0;
}
