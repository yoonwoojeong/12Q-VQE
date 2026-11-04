/*
 * Local VQE:
 *  - Exact reverse-mode adjoint differentiation for Trotterised Pauli-rotation circuits
 *  - L-BFGS (NLOPT_LD_LBFGS) and COBYLA (NLOPT_LN_COBYLA) optimisation
 *  - ADAPT-VQE operator pool growth via exact commutator gradients
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <math.h>
#include <nlopt.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const UCCSD *u;
    const PauliSum *ham;
    StateVec *sv;
    StateVec *lam;
    int nfev;
    int verbose;
    double best;
} ObjData;

static const double complex IPOW4_G[4] = {1.0, I, -1.0, -I};

/* <lam | i P(x,z) | psi> where P(x,z) = i^k X^x Z^z is Hermitian;
   real part is -2 Im <lam | P(x,z) | psi> for real-expectation derivative of <psi|H|psi> */
static double grad_pauli_term(const StateVec *lam, const StateVec *psi, uint32_t x, uint32_t z)
{
    int k = __builtin_popcount(x & z) & 3;
    double complex ik = IPOW4_G[k];
    double s_im = 0.0;
    size_t dim = psi->dim;
    #pragma omp parallel for reduction(+:s_im) if(dim >= 16384)
    for (size_t i = 0; i < dim; i++) {
        size_t j = i ^ x;
        double sg = (__builtin_popcount((unsigned)(j & z)) & 1) ? -1.0 : 1.0;
        double complex val = conj(lam->a[i]) * ik * sg * psi->a[j];
        s_im += cimag(val);
    }
    /* d/dphi <psi(phi)|H|psi(phi)> = 2 Re <H psi | i P psi> = -2 Im <lam | P | psi> */
    return -2.0 * s_im;
}

double vqe_energy(const UCCSD *u, const PauliSum *ham, const double *theta, StateVec *scratch)
{
    uccsd_prepare(u, theta, scratch);
    return sv_expect_sum(scratch, ham);
}

/* Reverse-mode adjoint differentiation over the UCCSD excitation sequence */
double vqe_energy_grad(const UCCSD *u, const PauliSum *ham, const double *theta,
                       double *grad, StateVec *psi, StateVec *lam)
{
    uccsd_prepare(u, theta, psi);
    sv_apply_sum(psi, ham, lam);
    double complex e = 0.0;
    for (size_t i = 0; i < psi->dim; i++) e += conj(psi->a[i]) * lam->a[i];

    for (int k = u->nexc - 1; k >= 0; k--) {
        const Excitation *ex = &u->exc[k];
        double gk = 0.0;
        for (int t = ex->nterm - 1; t >= 0; t--) {
            uint32_t x = ex->terms[t].x, z = ex->terms[t].z;
            double d = ex->terms[t].d;
            double phi = theta[k] * d;
            gk += d * grad_pauli_term(lam, psi, x, z);
            sv_apply_pauli_rotation(psi, x, z, -phi);
            sv_apply_pauli_rotation(lam, x, z, -phi);
        }
        grad[k] = gk;
    }
    return creal(e);
}

static double objective(unsigned n, const double *x, double *grad, void *data)
{
    (void)n;
    ObjData *d = data;
    double e;
    if (grad) {
        e = vqe_energy_grad(d->u, d->ham, x, grad, d->sv, d->lam);
    } else {
        e = vqe_energy(d->u, d->ham, x, d->sv);
    }
    d->nfev++;
    if (e < d->best) d->best = e;
    if (d->verbose && (d->nfev % (grad ? 10 : 50) == 0))
        printf("  VQE eval %5d  E = %.12f  (best %.12f)\n", d->nfev, e, d->best);
    return e;
}

int vqe_optimize_local(const UCCSD *u, const PauliSum *ham, const VQEOptions *opt, VQEResult *r)
{
    memset(r, 0, sizeof *r);
    r->nparam = u->nexc;
    r->theta = qc_calloc((size_t)u->nexc, sizeof(double));

    StateVec sv, lam;
    sv_init(&sv, u->nq);
    sv_init(&lam, u->nq);
    ObjData d = {u, ham, &sv, &lam, 0, opt->verbose, 1e300};

    nlopt_algorithm alg = opt->use_cobyla ? NLOPT_LN_COBYLA : NLOPT_LD_LBFGS;
    nlopt_opt o = nlopt_create(alg, (unsigned)u->nexc);
    if (!o) qc_die("nlopt_create failed");
    nlopt_set_min_objective(o, objective, &d);
    nlopt_set_xtol_rel(o, opt->xtol_rel > 0 ? opt->xtol_rel : 1e-10);
    nlopt_set_ftol_abs(o, 1e-13);
    nlopt_set_maxeval(o, opt->maxeval);
    if (opt->use_cobyla) nlopt_set_initial_step1(o, opt->init_step > 0 ? opt->init_step : 0.1);

    double fmin = 0.0;
    double t0 = qc_wall_time();
    int rc = nlopt_optimize(o, r->theta, &fmin);
    r->seconds = qc_wall_time() - t0;
    nlopt_destroy(o);

    r->energy = vqe_energy(u, ham, r->theta, &sv);   /* re-evaluate at the returned point */
    r->nfev = d.nfev;
    r->converged = (rc == NLOPT_XTOL_REACHED || rc == NLOPT_FTOL_REACHED || rc == NLOPT_SUCCESS);
    if (rc < 0) qc_warn("NLopt returned %d", rc);
    sv_free(&sv);
    sv_free(&lam);
    return rc;
}

/* ---- ADAPT-VQE (operator-by-operator ansatz construction with L-BFGS) ---- */

typedef struct {
    const UCCSD *pool;
    const int *op_idx;
    int depth;
    const PauliSum *ham;
    StateVec *psi;
    StateVec *lam;
    int nfev;
} AdaptObj;

static void adapt_prepare(const UCCSD *pool, const int *op_idx, int depth, const double *theta, StateVec *psi)
{
    sv_set_basis(psi, pool->hf_index);
    for (int m = 0; m < depth; m++) {
        const Excitation *ex = &pool->exc[op_idx[m]];
        for (int t = 0; t < ex->nterm; t++)
            sv_apply_pauli_rotation(psi, ex->terms[t].x, ex->terms[t].z, theta[m] * ex->terms[t].d);
    }
}

static double adapt_objective(unsigned n, const double *x, double *grad, void *data)
{
    AdaptObj *ao = data;
    ao->nfev++;
    adapt_prepare(ao->pool, ao->op_idx, (int)n, x, ao->psi);
    if (!grad) return sv_expect_sum(ao->psi, ao->ham);

    sv_apply_sum(ao->psi, ao->ham, ao->lam);
    double complex e = 0.0;
    for (size_t i = 0; i < ao->psi->dim; i++) e += conj(ao->psi->a[i]) * ao->lam->a[i];

    for (int m = (int)n - 1; m >= 0; m--) {
        const Excitation *ex = &ao->pool->exc[ao->op_idx[m]];
        double gm = 0.0;
        for (int t = ex->nterm - 1; t >= 0; t--) {
            uint32_t px = ex->terms[t].x, pz = ex->terms[t].z;
            double d = ex->terms[t].d;
            double phi = x[m] * d;
            gm += d * grad_pauli_term(ao->lam, ao->psi, px, pz);
            sv_apply_pauli_rotation(ao->psi, px, pz, -phi);
            sv_apply_pauli_rotation(ao->lam, px, pz, -phi);
        }
        grad[m] = gm;
    }
    return creal(e);
}

int adapt_vqe_optimize(const UCCSD *pool, const PauliSum *ham, int max_ops, double grad_tol,
                       int verbose, AdaptResult *res, Circuit *circ_out)
{
    memset(res, 0, sizeof *res);
    if (max_ops <= 0) max_ops = pool->nexc;
    if (grad_tol <= 0) grad_tol = 1e-4;
    res->op_idx = qc_calloc((size_t)max_ops, sizeof(int));
    res->theta = qc_calloc((size_t)max_ops, sizeof(double));

    StateVec psi, lam;
    sv_init(&psi, pool->nq);
    sv_init(&lam, pool->nq);
    double t0 = qc_wall_time();

    sv_set_basis(&psi, pool->hf_index);
    res->energy = sv_expect_sum(&psi, ham);

    for (int step = 0; step < max_ops; step++) {
        adapt_prepare(pool, res->op_idx, res->n_ops, res->theta, &psi);
        sv_apply_sum(&psi, ham, &lam);

        /* Pool commutator gradient: g_k = d/dtheta_k <psi|e^{-theta_k G_k} H e^{theta_k G_k}|psi>|_{0} */
        int best_k = -1;
        double max_g = 0.0, g_norm2 = 0.0;
        for (int k = 0; k < pool->nexc; k++) {
            const Excitation *ex = &pool->exc[k];
            double gk = 0.0;
            for (int t = 0; t < ex->nterm; t++)
                gk += ex->terms[t].d * grad_pauli_term(&lam, &psi, ex->terms[t].x, ex->terms[t].z);
            g_norm2 += gk * gk;
            if (fabs(gk) > max_g) {
                max_g = fabs(gk);
                best_k = k;
            }
        }
        res->max_grad = max_g;
        if (verbose)
            printf("  ADAPT step %2d  E = %.12f  max|g| = %.3e  ||g|| = %.3e  (op %d)\n",
                   step, res->energy, max_g, sqrt(g_norm2), best_k);
        if (max_g < grad_tol || best_k < 0) {
            res->converged = 1;
            break;
        }
        res->op_idx[res->n_ops] = best_k;
        res->theta[res->n_ops] = 0.0;
        res->n_ops++;

        AdaptObj ao = {pool, res->op_idx, res->n_ops, ham, &psi, &lam, 0};
        nlopt_opt o = nlopt_create(NLOPT_LD_LBFGS, (unsigned)res->n_ops);
        nlopt_set_min_objective(o, adapt_objective, &ao);
        nlopt_set_xtol_rel(o, 1e-10);
        nlopt_set_ftol_abs(o, 1e-13);
        nlopt_set_maxeval(o, 500);
        double fmin = res->energy;
        nlopt_optimize(o, res->theta, &fmin);
        nlopt_destroy(o);
        res->nfev += ao.nfev;
        adapt_prepare(pool, res->op_idx, res->n_ops, res->theta, &psi);
        res->energy = sv_expect_sum(&psi, ham);
    }
    res->seconds = qc_wall_time() - t0;

    if (circ_out) {
        circ_init(circ_out, pool->nq);
        for (int q = 0; q < pool->nq; q++)
            if ((pool->hf_mask >> q) & 1) circ_add(circ_out, G_X, q, -1, 0);
        for (int m = 0; m < res->n_ops; m++) {
            const Excitation *ex = &pool->exc[res->op_idx[m]];
            for (int t = 0; t < ex->nterm; t++)
                circ_append_pauli_rotation(circ_out, ex->terms[t].x, ex->terms[t].z,
                                           res->theta[m] * ex->terms[t].d, pool->nq);
        }
    }
    sv_free(&psi);
    sv_free(&lam);
    return 0;
}

void adapt_result_free(AdaptResult *r)
{
    free(r->op_idx);
    free(r->theta);
    memset(r, 0, sizeof *r);
}

void vqe_result_free(VQEResult *r)
{
    free(r->theta);
    memset(r, 0, sizeof *r);
}
