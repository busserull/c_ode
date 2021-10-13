#include "integrator.h"
#include <stdlib.h>
#include <stdarg.h>

static void weighted_sum(
    Vector * p_out,
    double time_scale,
    const Vector * p_x,
    const double * p_ab,
    const Vector * p_k,
    int number_of_ks
){
    for(int dim = 0; dim < p_x->dim; dim++){
        p_out->data[dim] = 0.0;

        for(int i = 0; i < number_of_ks; i++){
            p_out->data[dim] += p_ab[i] * p_k->data[i];
        }
        p_out->data[dim] *= time_scale;

        p_out->data[dim] += p_x->data[dim];
    }
}


static void table_allocate_space(IntegratorTable * p_table, int steps){
    /* Space requirements (n = steps) */
    /* a: n * (n - 1) / 2 */
    /* b: n */
    /* c: n */
    /* e: n */
    /* "implicit" zero: 1 */

    int size_a = (steps * steps - steps) / 2 + 1;
    int size_total = size_a + 3 * steps;
    double * p_base = (double *)malloc(size_total * sizeof(double));

    p_table->butcher_a = p_base;
    p_table->butcher_b = p_base + size_a + 1;
    p_table->butcher_c = p_base + size_a + steps + 1;
    p_table->butcher_e = p_base + size_a + 2 * steps + 1;

    p_table->steps = steps;
}

static void table_free_space(IntegratorTable * p_table){
    free(p_table->butcher_a);

    p_table->butcher_a = NULL;
    p_table->butcher_b = NULL;
    p_table->butcher_c = NULL;
    p_table->butcher_e = NULL;

    p_table->steps = 0;
}

static void table_set_a(IntegratorTable * p_table, ...){
    va_list ap;
    va_start(ap, p_table);

    int steps = p_table->steps;
    int size_a = (steps * steps - steps) / 2 + 1;

    p_table->butcher_a[0] = 0.0;
    for(int i = 1; i < size_a; i++){
        p_table->butcher_a[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void table_set_b(IntegratorTable * p_table, ...){
    va_list ap;
    va_start(ap, p_table);

    for(int i = 0; i < p_table->steps; i++){
        p_table->butcher_b[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void table_set_c(IntegratorTable * p_table, ...){
    va_list ap;
    va_start(ap, p_table);

    for(int i = 0; i < p_table->steps; i++){
        p_table->butcher_c[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void table_set_e(IntegratorTable * p_table, ...){
    va_list ap;
    va_start(ap, p_table);

    for(int i = 0; i < p_table->steps; i++){
        p_table->butcher_e[i] = va_arg(ap, double);
    }

    va_end(ap);
}

static void table_unset_e(IntegratorTable * p_table){
    p_table->butcher_e = NULL;
}

static IntegratorTable table_new_rk4(){
    IntegratorTable table;

    table_allocate_space(&table, 4);
    table_set_a(
        &table,
        1.0/2,
        0.0, 1.0/2,
        0.0, 0.0, 1.0
    );
    table_set_b(&table, 1.0/6, 1.0/3, 1.0/3, 1.0/6);
    table_set_c(&table, 0.0, 1.0/2, 1.0/2, 1.0);
    table_unset_e(&table);

    return table;
}


static Vector * scratchpad_new(int plant_dimension, int integrator_steps){
    Vector * s = (Vector *)malloc((integrator_steps + 1) * sizeof(Vector));

    for(int i = 0; i < integrator_steps + 1; i++){
        s[i].data = malloc(plant_dimension * sizeof(double));
        s[i].dim = plant_dimension;
    }

    return s;
}

static void scratchpad_delete(Vector * scratchpad, int integrator_steps){
    for(int i = 0; i < integrator_steps + 1; i++){
        free(scratchpad[i].data);
        scratchpad[i].dim = 0;
    }

    free(scratchpad);
}


Integrator integrator_new(Plant * p_plant, IntegratorMethod method){
    Integrator s;
    s.p_plant = p_plant;

    switch(method){
        case INTEGRATOR_METHOD_EXPLICIT_EULER:
            break;

        case INTEGRATOR_METHOD_RK4:
            s.table = table_new_rk4();
            break;

        case INTEGRATOR_METHOD_DORMAND_PRINCE:
            break;

        default:
            break;
    }

    s.scratchpad = scratchpad_new(p_plant->dim, s.table.steps);

    return s;
}

void integrator_step(
    Integrator * p_integrator,
    Vector * p_xout,
    double step_size,
    double t,
    const Vector * p_x,
    const Vector * p_u
){
    Vector * p_xtemp = p_integrator->scratchpad;
    Vector * p_k = p_integrator->scratchpad + 1;

    int a_offset = 0;
    for(int i = 0; i < p_integrator->table.steps; i++){
        weighted_sum(
            p_xtemp,
            step_size,
            p_x,
            p_integrator->table.butcher_a + a_offset + 1,
            p_k,
            i
        );

        p_integrator->p_plant->xdot(
            p_integrator->p_plant->p_params,
            p_k + i,
            t + p_integrator->table.butcher_c[i] * step_size,
            p_xtemp,
            p_u
        );

        a_offset += i;
    }

    weighted_sum(
        p_xtemp,
        step_size,
        p_x,
        p_integrator->table.butcher_b,
        p_k,
        p_integrator->table.steps - 1
    );

    for(int i = 0; i < p_x->dim; i++){
        p_xout->data[i] = p_x->data[i];
    }
}

void integrator_delete(Integrator * p_integrator){
    scratchpad_delete(p_integrator->scratchpad, p_integrator->table.steps);
    table_free_space(&(p_integrator->table));
}




/* Solver solver_new(Plant * p_plant){ */
/*     Solver s; */

/*     s.p_plant = p_plant; */
/*     s.scratchpad = (Vector *)malloc(8 * sizeof(Vector)); */

/*     for(int i = 0; i < p_plant->dim + 1; i++){ */
/*         s.scratchpad[i].dim = p_plant->dim; */
/*         s.scratchpad[i].data = malloc(p_plant->dim * sizeof(double)); */
/*     } */

/*     return s; */
/* } */

/* void solver_delete(Solver * p_solver){ */
/*     for(int i = 0; i < p_solver->p_plant->dim + 1; i++){ */
/*         free(p_solver->scratchpad[i].data); */
/*         p_solver->scratchpad[i].dim = 0; */
/*     } */

/*     free(p_solver->scratchpad); */
/*     p_solver->scratchpad = NULL; */
/* } */

/* void solver_step( */
/*     Solver * p_solver, */
/*     Vector * p_xout, */
/*     double t, */
/*     double dt, */
/*     const Vector * p_x, */
/*     const Vector * p_u */
/* ){ */
/*     static const double butcher_a[7][6] = { */
/*         {0.0}, */
/*         {1.0/5}, */
/*         {3.0/40, 9.0/40}, */
/*         {44.0/45, -56.0/15, 32.0/9}, */
/*         {19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729}, */
/*         {9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656}, */
/*         {35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84} */
/*     }; */

/*     static const double butcher_b[7] = { */
/*         35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0.0 */
/*     }; */

/*     static const double butcher_c[7] = { */
/*         0.0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1.0, 1.0 */
/*     }; */

/*     int plant_dim = p_solver->p_plant->dim; */
/*     PlantEq f = p_solver->p_plant->xdot; */

/*     Vector * p_xtemp = p_solver->scratchpad; */
/*     Vector * k = p_solver->scratchpad + 1; */

/*     for(int i = 0; i < 7; i++){ */
/*         sum_x_k(p_xtemp, p_x, k, butcher_a[i], i); */

/*         f( */
/*             k + i, */
/*             t + dt * butcher_c[i], */
/*             p_xtemp, */
/*             p_u */
/*         ); */
/*     } */

/*     for(int dim = 0; dim < plant_dim; dim++){ */
/*         p_xtemp->data[dim] = 0.0; */

/*         for(int i = 0; i < 7; i++){ */
/*             p_xtemp->data */
/*         } */
/*     } */
/* } */
