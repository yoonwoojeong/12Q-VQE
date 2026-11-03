/*
 * Determinant-based full CI in the active space (classical reference).
 *
 * The Hamiltonian matrix is built directly in the Fock space from the
 * spin-orbital integrals using explicit creation/annihilation sign rules.
 * This does not use the Jordan-Wigner mapping, so comparing it with the
 * Pauli Hamiltonian is an independent check of the qubit operator.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <lapacke.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* apply a_q (ann=1) or a+_q (ann=0) to determinant *det, returns sign, 0 if zero */
static int apply_op(uint32_t *det, int q, int ann)
{
    uint32_t bit = 1u << q;
    int occ = (*det & bit) != 0;
    if (ann && !occ) return 0;
    if (!ann && occ) return 0;
    int sign = (__builtin_popcount(*det & (bit - 1)) & 1) ? -1 : 1;
    *det ^= bit;
    return sign;
}

static void enumerate_dets(int ns, int nel, uint32_t **out, int *nout)
{
    int cap = 1;
    /* binomial upper bound */
    long long c = 1;
    for (int k = 0; k < nel; k++) c = c * (ns - k) / (k + 1);
    cap = (int)c;
    uint32_t *d = qc_calloc((size_t)cap, sizeof(uint32_t));
    int n = 0;
    for (uint32_t m = 0; m < (1u << ns); m++)
        if (__builtin_popcount(m) == nel) d[n++] = m;
    *out = d;
    *nout = n;
}

static int find_det(const uint32_t *dets, int n, uint32_t m)
{
    for (int i = 0; i < n; i++)
        if (dets[i] == m) return i;
    return -1;
}

int fci_solve(FCIResult *r, const FermionHam *fh, const ActiveSpace *as, int nelec)
{
    memset(r, 0, sizeof *r);
    int ns = fh->nspin;
    int nd;
    enumerate_dets(ns, nelec, &r->dets, &nd);
    r->ndet = nd;
    double *H = qc_calloc((size_t)nd * nd, sizeof(double));

    for (int J = 0; J < nd; J++) {
        /* one-body: sum h_pq a+_p a_q |J> */
        for (int q = 0; q < ns; q++)
            for (int p = 0; p < ns; p++) {
                double h = fh->h_spin[(size_t)p * ns + q];
                if (h == 0.0) continue;
                uint32_t d = r->dets[J];
                int s1 = apply_op(&d, q, 1);
                if (!s1) continue;
                int s2 = apply_op(&d, p, 0);
                if (!s2) continue;
                int ii = find_det(r->dets, nd, d);
                if (ii >= 0) H[(size_t)ii * nd + J] += (double)(s1 * s2) * h;
            }
        /* two-body: 1/2 sum <pq|rs> a+_p a+_q a_s a_r |J> */
        for (int p = 0; p < ns; p++)
            for (int q = 0; q < ns; q++)
                for (int rr = 0; rr < ns; rr++)
                    for (int s = 0; s < ns; s++) {
                        double g = fh->g_spin[(((size_t)p * ns + q) * ns + rr) * ns + s];
                        if (g == 0.0) continue;
                        uint32_t d = r->dets[J];
                        int s1 = apply_op(&d, rr, 1);
                        if (!s1) continue;
                        int s2 = apply_op(&d, s, 1);
                        if (!s2) continue;
                        int s3 = apply_op(&d, q, 0);
                        if (!s3) continue;
                        int s4 = apply_op(&d, p, 0);
                        if (!s4) continue;
                        int ii = find_det(r->dets, nd, d);
                        if (ii >= 0) H[(size_t)ii * nd + J] += 0.5 * (double)(s1 * s2 * s3 * s4) * g;
                    }
        H[(size_t)J * nd + J] += as->ecore;
    }
    /* symmetry check */
    double asym = 0.0;
    for (int i = 0; i < nd; i++)
        for (int j = 0; j < nd; j++) {
            double d = fabs(H[(size_t)i * nd + j] - H[(size_t)j * nd + i]);
            if (d > asym) asym = d;
        }
    if (asym > 1e-10) qc_warn("FCI matrix asymmetry %.3e", asym);

    double *w = qc_calloc((size_t)nd, sizeof(double));
    if (LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nd, H, nd, w) != 0) qc_die("dsyev failed in FCI");
    r->eigvals = w;
    r->ground = w[0];
    r->ground_vec = qc_calloc((size_t)nd, sizeof(double));
    for (int i = 0; i < nd; i++) r->ground_vec[i] = H[(size_t)i * nd + 0];
    free(H);
    return 0;
}

double fci_pauli_restricted_ground(const FermionHam *fh, const FCIResult *fci)
{
    int ns = fh->nspin;
    int nd = fci->ndet;
    double *M = qc_calloc((size_t)nd * nd, sizeof(double));
    StateVec in, out;
    sv_init(&in, ns);
    sv_init(&out, ns);
    for (int J = 0; J < nd; J++) {
        sv_set_basis(&in, fci->dets[J]);
        sv_apply_sum(&in, &fh->ham, &out);
        for (int ii = 0; ii < nd; ii++) {
            double complex a = out.a[fci->dets[ii]];
            M[(size_t)ii * nd + J] = creal(a);
            if (fabs(cimag(a)) > 1e-10) qc_warn("Pauli-restricted FCI: imaginary matrix element");
        }
        /* leakage outside the N-electron sector must vanish */
        double leak = 0.0;
        for (size_t k = 0; k < out.dim; k++) {
            int in_sector = 0;
            for (int ii = 0; ii < nd; ii++)
                if (fci->dets[ii] == k) { in_sector = 1; break; }
            if (!in_sector) leak += cabs(out.a[k]);
        }
        if (leak > 1e-9) qc_warn("Pauli Hamiltonian leaks out of the particle-number sector (%.3e)", leak);
    }
    double *w = qc_calloc((size_t)nd, sizeof(double));
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', nd, M, nd, w);
    double g = w[0];
    free(w);
    free(M);
    sv_free(&in);
    sv_free(&out);
    return g;
}

void fci_free(FCIResult *r)
{
    free(r->eigvals);
    free(r->ground_vec);
    free(r->dets);
    memset(r, 0, sizeof *r);
}
