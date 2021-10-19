#include "types.h"
#include "rkdp.h"

#include <stdio.h>
#include <stdlib.h>

void mass_damper_spring(
    void * p_plant_params,
    Vector x_dot,
    double t,
    const Vector x,
    const Vector u
){
    double m = ((double *)p_plant_params)[0];
    double d = ((double *)p_plant_params)[1];
    double k = ((double *)p_plant_params)[2];

    x_dot[0] = x[1];
    x_dot[1] = (-k * x[0] - d * x[1] + u[1]) / m;
}



int main(){
    Plant mds;
    double mds_params[3] = {1.0, 0.1, 1.0};
    mds.p_plant_params = (void *)mds_params;
    mds.plant_function = mass_damper_spring;

    void * p_rkdp_wa = rkdp_working_area_new(2);

    double x[2] = {1.0, 0.0};
    double u[2] = {0.0, 0.0};

    double t = 0.0;
    double h = 0.1;

    printf("x = [");
    while(t < 30.0){
        printf("%f,...\n", x[0]);
        rkdp_step(p_rkdp_wa, &mds, x, NULL, t, x, u, h);
        t += h;
    }
    printf("]; plot(x);\n");

    rkdp_working_area_delete(p_rkdp_wa);
    return 0;
}
