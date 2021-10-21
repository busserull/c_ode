#include "types.h"
#include "rkdp.h"

static RKDP_WA_NEW(m_rkdp_wa, 2);

static double x[2] = {1.0, 0.0};
static double u[2] = {0.0, 0.0};

static double t = 0.0;
static double h = 0.1;

static double plant_params[3] = {1.0, 0.1, 1.0};

static void plant_function(
    double * p_plant_params,
    Vector x_dot,
    double t,
    const Vector x,
    const Vector u
){
    double m = p_plant_params[0];
    double d = p_plant_params[1];
    double k = p_plant_params[2];

    x_dot[0] = x[1];
    x_dot[1] = (-k * x[0] - d * x[1] + u[1]) / m;
}

static Plant plant = {
    plant_params,
    plant_function
};


double get_t(){
    return t;
}

double get_x(int index){
    return x[index];
}

double get_u(int index){
    return u[index];
}

double get_plant_param(int index){
    return plant_params[index];
}


void set_t(double new_t){
    t = new_t;
}

void set_x(int index, double new_x){
    x[index] = new_x;
}

void set_u(int index, double new_u){
    u[index] = new_u;
}

void set_plant_param(int index, double new_param){
    plant_params[index] = new_param;
}


double step(){
    rkdp_step(m_rkdp_wa, &plant, x, NULL, t, x, u, h);
    t = t + h;

    return t;
}
