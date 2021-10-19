#ifndef TYPES_H
#define TYPES_H

typedef double * Vector;

typedef void (* PlantFunction)(
    void * p_plant_params,
    Vector x_dot,
    double t,
    const Vector x,
    const Vector u
);

typedef struct {
    void * p_plant_params;
    PlantFunction plant_function;
} Plant;

#endif
