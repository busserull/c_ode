/* #include "integrator.h" */
#include "solver.h"
#include <stdlib.h>
#include <stdarg.h>

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

typedef struct {
    int steps;
    double * butcher_a;
    double * butcher_b;
    double * butcher_c;
    double * butcher_e;
} IntegratorMethod;

static void method_allocate_space(IntegratorMethod * p_method, int steps){
    /* Space requirements (n = steps) */
    /* a: n * (n - 1) / 2 */
    /* b: n */
    /* c: n */
    /* e: n */
    /* "implicit" zero: 1 */

    int size_a = (steps * steps - steps) / 2 + 1;
    int size_total = size_a + 3 * steps;
    double * p_base = (double *)malloc(size_total * sizeof(double));

    p_method->butcher_a = p_base;
    p_method->butcher_b = p_base + size_a + 1;
    p_method->butcher_c = p_base + size_a + steps + 1;
    p_method->butcher_e = p_base + size_a + 2 * steps + 1;

    p_method->steps = steps;
}

static void method_free_space(IntegratorMethod * p_method){
    free(p_method->butcher_a);

    p_method->butcher_a = NULL;
    p_method->butcher_b = NULL;
    p_method->butcher_c = NULL;
    p_method->butcher_e = NULL;

    p_method->steps = 0;
}

static void method_set_a(IntegratorMethod * p_method, ...){
    va_list ap;
    va_start(ap, p_method);

    int steps = p_method->steps;
    int size_a = (steps * steps - steps) / 2 + 1;

    p_method->butcher_a[0] = 0.0;
    for(int i = 1; i < size_a; i++){
        p_method->butcher_a[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void method_set_b(IntegratorMethod * p_method, ...){
    va_list ap;
    va_start(ap, p_method);

    for(int i = 0; i < p_method->steps; i++){
        p_method->butcher_b[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void method_set_c(IntegratorMethod * p_method, ...){
    va_list ap;
    va_start(ap, p_method);

    for(int i = 0; i < p_method->steps; i++){
        p_method->butcher_c[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void method_set_e(IntegratorMethod * p_method, ...){
    va_list ap;
    va_start(ap, p_method);

    for(int i = 0; i < p_method->steps; i++){
        p_method->butcher_e[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void method_unset_e(IntegratorMethod * p_method){
    p_method->butcher_e = NULL;
}

static IntegratorMethod method_new_rk4(){
    IntegratorMethod method;

    method_allocate_space(&method, 4);
    method_set_a(
        &method,
        1.0/2,
        0.0, 1.0/2,
        0.0, 0.0, 1.0
    );
    method_set_b(&method, 1.0/6, 1.0/3, 1.0/3, 1.0/6);
    method_set_c(&method, 0.0, 1.0/2, 1.0/2, 1.0);
    method_unset_e(&method);

    return method;
}

Integrator integrator_new(Plant * p_plant){
    Integrator s;

    return s;
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

    for(int dim = 0; dim < plant_dim; dim++){
        p_xtemp->data[dim] = 0.0;

        for(int i = 0; i < 7; i++){
            p_xtemp->data
        }
    }
}
