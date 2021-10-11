#ifndef SOLVER_H
#define SOLVER_H
#include "types.h"

typedef enum {
    SOLVER_METHOD_EXPLICIT_EULER,
    SOLVER_METHOD_RK4,
    SOLVER_METHOD_DORMAND_PRINCE
} SolverMethod;

typedef struct {
    Plant * p_plant;
    Vector * scratchpad;
    void * p_method;
} Solver;

Solver solver_new(Plant * p_plant);

void solver_delete(Solver * p_solver);

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

typedef struct {
    Plant * p_plant;
    Vector * scratchpad;
} Integrator;

Integrator integrator_new(Plant * p_plant);

#endif
