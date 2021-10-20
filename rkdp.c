#include "rkdp.h"
#include <stdlib.h>

static void vector_null(Vector v, int dim){
    for(int i = 0; i < dim; i++){
        v[i] = 0.0;
    }
}

static void vector_add_scale(Vector v, int dim, double scale, const Vector k){
    for(int i = 0; i < dim; i++){
        v[i] += scale * k[i];
    }
}

static void vector_scale(Vector v, int dim, double scale){
    for(int i = 0; i < dim; i++){
        v[i] *= scale;
    }
}

static void vector_add(Vector v, int dim, const Vector x){
    for(int i = 0; i < dim; i++){
        v[i] += x[i];
    }
}

void rkdp_step(
    double * p_working_area,
    Plant * p_plant,
    Vector x_out,
    Vector e_out,
    double t,
    const Vector x,
    const Vector u,
    double step_size
){
    int dim = *((int *)p_working_area);

    double * x_temp = p_working_area + 1;
    double * k0 = x_temp + dim;
    double * k1 = x_temp + 2 * dim;
    double * k2 = x_temp + 3 * dim;
    double * k3 = x_temp + 4 * dim;
    double * k4 = x_temp + 5 * dim;
    double * k5 = x_temp + 6 * dim;
    double * k6 = x_temp + 7 * dim;

    void * p_plant_params = p_plant->p_plant_params;
    PlantFunction f = p_plant->plant_function;

    double h = step_size;

    /* Step 0 */
    f(p_plant_params, k0, t, x, u);

    /* Step 1 */
    double c1 = 1.0/5.0;
    double a10 = 1.0/5.0;
    vector_null(x_temp, dim);
    vector_add_scale(x_temp, dim, a10, k0);
    vector_scale(x_temp, dim, h);
    vector_add(x_temp, dim, x);

    f(p_plant_params, k1, t + c1 * h, x_temp, u);

    /* Step 2 */
    double c2 = 3.0/10.0;
    double a20 = 3.0/40.0;
    double a21 = 9.0/40.0;
    vector_null(x_temp, dim);
    vector_add_scale(x_temp, dim, a20, k0);
    vector_add_scale(x_temp, dim, a21, k1);
    vector_scale(x_temp, dim, h);
    vector_add(x_temp, dim, x);

    f(p_plant_params, k2, t + c2 * h, x_temp, u);

    /* Step 3 */
    double c3 = 4.0/5.0;
    double a30 = 44.0/45.0;
    double a31 = -56.0/15.0;
    double a32 = 32.0/9.0;
    vector_null(x_temp, dim);
    vector_add_scale(x_temp, dim, a30, k0);
    vector_add_scale(x_temp, dim, a31, k1);
    vector_add_scale(x_temp, dim, a32, k2);
    vector_scale(x_temp, dim, h);
    vector_add(x_temp, dim, x);

    f(p_plant_params, k3, t + c3 * h, x_temp, u);

    /* Step 4 */
    double c4 = 8.0/9.0;
    double a40 = 19372.0/6561.0;
    double a41 = -25360.0/2187.0;
    double a42 = 64448.0/6561.0;
    double a43 = -212.0/729.0;
    vector_null(x_temp, dim);
    vector_add_scale(x_temp, dim, a40, k0);
    vector_add_scale(x_temp, dim, a41, k1);
    vector_add_scale(x_temp, dim, a42, k2);
    vector_add_scale(x_temp, dim, a43, k3);
    vector_scale(x_temp, dim, h);
    vector_add(x_temp, dim, x);

    f(p_plant_params, k4, t + c4 * h, x_temp, u);

    /* Step 5 */
    double c5 = 1.0;
    double a50 = 9017.0/3168.0;
    double a51 = -355.0/33.0;
    double a52 = 46732.0/5247.0;
    double a53 = 49.0/176.0;
    double a54 = -5103.0/18656.0;
    vector_null(x_temp, dim);
    vector_add_scale(x_temp, dim, a50, k0);
    vector_add_scale(x_temp, dim, a51, k1);
    vector_add_scale(x_temp, dim, a52, k2);
    vector_add_scale(x_temp, dim, a53, k3);
    vector_add_scale(x_temp, dim, a54, k4);
    vector_scale(x_temp, dim, h);
    vector_add(x_temp, dim, x);

    f(p_plant_params, k5, t + c5 * h, x_temp, u);

    /* Step 6 */
    if(e_out != NULL){
        double c6 = 1.0;
        double a60 = 35.0/384.0;
        double a62 = 500.0/1113.0;
        double a63 = 125.0/192.0;
        double a64 = -2187.0/6784.0;
        double a65 = 11.0/84.0;
        vector_null(x_temp, dim);
        vector_add_scale(x_temp, dim, a60, k0);
        vector_add_scale(x_temp, dim, a62, k2);
        vector_add_scale(x_temp, dim, a63, k3);
        vector_add_scale(x_temp, dim, a64, k4);
        vector_add_scale(x_temp, dim, a65, k5);
        vector_scale(x_temp, dim, h);
        vector_add(x_temp, dim, x);

        f(p_plant_params, k6, t + c6 * h, x_temp, u);
    }

    /* Output */
    double b0 = 35.0/384.0;
    double b2 = 500.0/1113.0;
    double b3 = 125.0/192.0;
    double b4 = -2187.0/6784.0;
    double b5 = 11.0/84.0;
    vector_null(x_temp, dim);
    vector_add_scale(x_temp, dim, b0, k0);
    vector_add_scale(x_temp, dim, b2, k2);
    vector_add_scale(x_temp, dim, b3, k3);
    vector_add_scale(x_temp, dim, b4, k4);
    vector_add_scale(x_temp, dim, b5, k5);
    vector_scale(x_temp, dim, h);
    vector_add(x_temp, dim, x);

    /* Error estimate */
    if(e_out != NULL){
        double e0 = 5179.0/57600.0;
        double e2 = 7571.0/16695.0;
        double e3 = 393.0/640.0;
        double e4 = -92097.0/339200.0;
        double e5 = 187.0/2100.0;
        double e6 = 1.0/40.0;
        vector_null(e_out, dim);
        vector_add_scale(e_out, dim, e0, k0);
        vector_add_scale(e_out, dim, e2, k2);
        vector_add_scale(e_out, dim, e3, k3);
        vector_add_scale(e_out, dim, e4, k4);
        vector_add_scale(e_out, dim, e5, k5);
        vector_add_scale(e_out, dim, e6, k6);
        vector_scale(e_out, dim, h);
        vector_add(e_out, dim, x);

        for(int i = 0; i < dim; i++){
            e_out[i] = x_temp[i] - e_out[i];
        }
    }

    /* Copy to destination last, in case x_out == x */
    for(int i = 0; i < dim; i++){
        x_out[i] = x_temp[i];
    }
}
