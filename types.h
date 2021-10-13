#ifndef TYPES_H
#define TYPES_H

typedef struct {
    int dim;
    double * data;
} Vector;

typedef void (* PlantEq)(
    void * p_plant_params,
    Vector * p_xdot,
    double t,
    const Vector * p_x,
    const Vector * p_u
);

typedef struct {
    int dim;
    void * p_params;
    PlantEq xdot;
} Plant;

#endif
