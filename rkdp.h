#ifndef RKDP_H
#define RKDP_H
#include "types.h"

#define RKDP_WA_NEW(name, plant_dim) \
    RKDP_WA (name)[1 + 8 * (plant_dim)] = \
    {{ .as_int = (plant_dim) }}

typedef union {
    int as_int;
    double as_double;
} RKDP_WA;

void rkdp_step(
    RKDP_WA * p_working_area,
    Plant * p_plant,
    Vector x_out,
    Vector e_out,
    double t,
    const Vector x,
    const Vector u,
    double step_size
);

#endif
