#include "solver.h"
#include <stdlib.h>

/* static scalar_product(Vector * out, ){ */

/* } */
static void sum_x_k(
    Vector * p_out,
    const Vector * p_x,
    const Vector * p_k,
    const double * a,
    int number_of_ks
){
    for(int i = 0; i < p_x->dim; i++){
        p_out->data[i] = p_x->data[i];
    }

    for(int k = 0; k < number_of_ks; k++){
        for(int i = 0; i < p_x->dim; i++){
            p_out->data[i] += a[i] * p_k[k].data[i];
        }
    }
}

Solver solver_new(Plant * p_plant){
    Solver s;

    s.p_plant = p_plant;
    s.scratchpad = (Vector *)malloc(8 * sizeof(Vector));

    for(int i = 0; i < p_plant->dim + 1; i++){
        s.scratchpad[i].dim = p_plant->dim;
        s.scratchpad[i].data = malloc(p_plant->dim * sizeof(double));
    }

    return s;
}

void solver_delete(Solver * p_solver){
    for(int i = 0; i < p_solver->p_plant->dim + 1; i++){
        free(p_solver->scratchpad[i].data);
        p_solver->scratchpad[i].dim = 0;
    }

    free(p_solver->scratchpad);
    p_solver->scratchpad = NULL;
}

void solver_step(
    Solver * p_solver,
    Vector * p_xout,
    double t,
    double dt,
    const Vector * p_x,
    const Vector * p_u
){
    static const double butcher_a[7][6] = {
        {0.0},
        {1.0/5},
        {3.0/40, 9.0/40},
        {44.0/45, -56.0/15, 32.0/9},
        {19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729},
        {9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656},
        {35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84}
    };

    static const double butcher_b[7] = {
        35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0.0
    };

    static const double butcher_c[7] = {
        0.0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1.0, 1.0
    };

    int plant_dim = p_solver->p_plant->dim;
    PlantEq f = p_solver->p_plant->xdot;

    Vector * p_xtemp = p_solver->scratchpad;
    Vector * k = p_solver->scratchpad + 1;

    for(int i = 0; i < 7; i++){
        sum_x_k(p_xtemp, p_x, k, butcher_a[i], i);

        f(
            k + i,
            t + dt * butcher_c[i],
            p_xtemp,
            p_u
        );
    }

}
