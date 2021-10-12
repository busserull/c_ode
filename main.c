#include "integrator.h"

#include <stdio.h>

/* void mass_damper_spring(Vector * p_xdot, const Vector * p_x){ */

/* } */

int main(){
    Plant p;
    p.dim = 2;
    p.p_params = NULL;
    p.xdot = NULL;

    Integrator s = integrator_new(&p, INTEGRATOR_METHOD_RK4);

    integrator_delete(&s);

    return 0;
}
