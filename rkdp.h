#ifndef RKDP_H
#define RKDP_H
#include "types.h"

void * rkdp_working_area_new(int plant_dimension);

void rkdp_working_area_delete(void * p_working_area);

void rkdp_step(
    void * p_rkdp_working_area,
    Plant * p_plant,
    Vector x_out,
    Vector e_out,
    double t,
    const Vector x,
    const Vector u,
    double step_size
);

#endif
