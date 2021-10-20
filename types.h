#ifndef TYPES_H
#define TYPES_H

typedef double * Vector;

typedef void (* PlantFunction)(
    double * p_plant_params,
    Vector x_dot,
    double t,
    const Vector x,
    const Vector u
);

typedef struct {
    double * p_plant_params;
    PlantFunction plant_function;
} Plant;

#endif
