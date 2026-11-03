/*
 * Local VQE: UCCSD on the 12-qubit statevector simulator, COBYLA (NLopt).
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <nlopt.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const UCCSD *u;
    const PauliSum *ham;
    StateVec *sv;
    int nfev;
    int verbose;
    double best;
} ObjData;

double vqe_energy(const UCCSD *u, const PauliSum *ham, const double *theta, StateVec *scratch)
{
    uccsd_prepare(u, theta, scratch);
    return sv_expect_sum(scratch, ham);
}

static double objective(unsigned n, const double *x, double *grad, void *data)
{
    (void)n;
    (void)grad;
    ObjData *d = data;
    double e = vqe_energy(d->u, d->ham, x, d->sv);
    d->nfev++;
    if (e < d->best) d->best = e;
    if (d->verbose && (d->nfev % 50 == 0))
        printf("  COBYLA eval %5d  E = %.12f  (best %.12f)\n", d->nfev, e, d->best);
    return e;
}

int vqe_optimize_local(const UCCSD *u, const PauliSum *ham, const VQEOptions *opt, VQEResult *r)
{
    memset(r, 0, sizeof *r);
    r->nparam = u->nexc;
    r->theta = qc_calloc((size_t)u->nexc, sizeof(double));

    StateVec sv;
    sv_init(&sv, u->nq);
    ObjData d = {u, ham, &sv, 0, opt->verbose, 1e300};

    nlopt_opt o = nlopt_create(NLOPT_LN_COBYLA, (unsigned)u->nexc);
    if (!o) qc_die("nlopt_create(COBYLA) failed");
    nlopt_set_min_objective(o, objective, &d);
    nlopt_set_xtol_rel(o, opt->xtol_rel);
    nlopt_set_maxeval(o, opt->maxeval);
    nlopt_set_initial_step1(o, opt->init_step);

    double fmin = 0.0;
    double t0 = qc_wall_time();
    int rc = nlopt_optimize(o, r->theta, &fmin);
    r->seconds = qc_wall_time() - t0;
    nlopt_destroy(o);

    r->energy = vqe_energy(u, ham, r->theta, &sv);   /* re-evaluate at the returned point */
    r->nfev = d.nfev;
    r->converged = (rc == NLOPT_XTOL_REACHED || rc == NLOPT_FTOL_REACHED || rc == NLOPT_SUCCESS);
    if (rc < 0) qc_warn("NLopt COBYLA returned %d", rc);
    sv_free(&sv);
    return rc;
}

void vqe_result_free(VQEResult *r)
{
    free(r->theta);
    memset(r, 0, sizeof *r);
}
