#ifndef RKDP_H
#define RKDP_H
#include "types.h"

#define RKDP_WA(name, plant_dim) \
    double (name)[1 + 8 * (plant_dim)]; \
    *((int *)(name)) = (plant_dim)

void rkdp_step(
    double * p_rkdp_working_area,
    Plant * p_plant,
    Vector x_out,
    Vector e_out,
    double t,
    const Vector x,
    const Vector u,
    double step_size
);

#endif
