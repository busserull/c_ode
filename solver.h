#ifndef SOLVER_H
#define SOLVER_H
#include "types.h"

typedef struct {
    Plant * p_plant;
    Vector * scratchpad;
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
