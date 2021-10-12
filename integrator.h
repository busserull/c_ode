#ifndef SOLVER_H
#define SOLVER_H
#include "types.h"

void solver_step(
    Solver * p_solver,
    Vector * p_xout,
    double t,
    double dt,
    const Vector * p_x,
    const Vector * p_u
);


#endif

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

typedef enum {
    INTEGRATOR_PARAMETERS_EXPLICIT_EULER,
    INTEGRATOR_PARAMETERS_RK4,
    INTEGRATOR_PARAMETERS_DORMAND_PRINCE,
} IntegratorMethod;

typedef struct {
    int steps;
    double * butcher_a;
    double * butcher_b;
    double * butcher_c;
    double * butcher_e;
} IntegratorTable;

typedef struct {
    Plant * p_plant;
    Vector * scratchpad;
    IntegratorTable table;
} Integrator;

Integrator integrator_new(Plant * p_plant, IntegratorMethod method);

void integrator_delete(Integrator * p_integrator);

#endif
