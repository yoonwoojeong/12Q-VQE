/*
 * Second-quantised active-space Hamiltonian -> Jordan-Wigner Pauli sum.
 *
 * Spin-orbital index s = 2p + sigma (sigma = 0 alpha, 1 beta), qubit s.
 *   H = E_core + sum_st h_st a+_s a_t + 1/2 sum_pqrs <pq|rs> a+_p a+_q a_s a_r
 *   <pq|rs> = (p_sp r_sr | q_sq s_ss) = (pr|qs) delta(sigma_p,sigma_r) delta(sigma_q,sigma_s)
 *   a_j  = (prod_{k<j} Z_k) (X_j + i Y_j) / 2
 *   a+_j = (prod_{k<j} Z_k) (X_j - i Y_j) / 2
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void fermion_ladder(int j, int dagger, PauliSum *out)
{
    /* Z string on modes < j */
    uint32_t zmask = (j > 0) ? ((1u << j) - 1u) : 0u;
    uint32_t xb = 1u << j;
    pauli_sum_init(out, j + 1);
    /* (X + i Y)/2 = 1/2 P(xb,0) + i/2 P(xb,xb) ;  (X - i Y)/2 = 1/2 P(xb,0) - i/2 P(xb,xb) */
    double complex sgn = dagger ? -I : I;
    pauli_sum_add(out, xb, 0, 0.5 + 0.0 * I);
    pauli_sum_add(out, xb, xb, 0.5 * sgn);
    /* prepend Z string: multiply by P(0, zmask) */
    if (zmask) {
        PauliSum zs;
        pauli_sum_init(&zs, j + 1);
        pauli_sum_add(&zs, 0, zmask, 1.0);
        PauliSum prod = pauli_mul(&zs, out, 1e-300);
        pauli_sum_free(&zs);
        pauli_sum_free(out);
        *out = prod;
    }
}

int fermion_build_hamiltonian(FermionHam *fh, const ActiveSpace *as)
{
    int norb = as->nact;
    int ns = 2 * norb;
    memset(fh, 0, sizeof *fh);
    fh->nspin = ns;
    fh->ecore = as->ecore;

    /* spin-orbital one- and two-body integrals */
    fh->h_spin = qc_calloc((size_t)ns * ns, sizeof(double));
    fh->g_spin = qc_calloc((size_t)ns * ns * ns * ns, sizeof(double));
    for (int s = 0; s < ns; s++)
        for (int t = 0; t < ns; t++) {
            if ((s & 1) != (t & 1)) continue;
            fh->h_spin[(size_t)s * ns + t] = as->h[(size_t)(s >> 1) * norb + (t >> 1)];
        }
    /* physicist <pq|rs> from chemist (pr|qs) */
    for (int p = 0; p < ns; p++)
        for (int q = 0; q < ns; q++)
            for (int r = 0; r < ns; r++)
                for (int s = 0; s < ns; s++) {
                    if ((p & 1) != (r & 1)) continue;
                    if ((q & 1) != (s & 1)) continue;
                    int P = p >> 1, R = r >> 1, Q = q >> 1, S = s >> 1;
                    double v = as->g[(((size_t)P * norb + R) * norb + Q) * norb + S];
                    fh->g_spin[(((size_t)p * ns + q) * ns + r) * ns + s] = v;
                }

    /* ladder operators */
    PauliSum *A = qc_calloc((size_t)ns, sizeof(PauliSum));
    PauliSum *C = qc_calloc((size_t)ns, sizeof(PauliSum));
    for (int j = 0; j < ns; j++) {
        fermion_ladder(j, 0, &A[j]);
        fermion_ladder(j, 1, &C[j]);
    }

    pauli_sum_init(&fh->ham, ns);
    pauli_sum_add(&fh->ham, 0, 0, fh->ecore + 0.0 * I);

    /* one-body */
    for (int s = 0; s < ns; s++)
        for (int t = 0; t < ns; t++) {
            double h = fh->h_spin[(size_t)s * ns + t];
            if (fabs(h) < 1e-14) continue;
            PauliSum prod = pauli_mul(&C[s], &A[t], 1e-15);
            pauli_sum_append(&fh->ham, &prod, h);
            pauli_sum_free(&prod);
        }

    /* two-body: 1/2 <pq|rs> C_p C_q A_s A_r */
    for (int p = 0; p < ns; p++)
        for (int q = 0; q < ns; q++) {
            if (p == q) continue;
            PauliSum CC = pauli_mul(&C[p], &C[q], 1e-15);
            for (int r = 0; r < ns; r++)
                for (int s = 0; s < ns; s++) {
                    if (r == s) continue;
                    double g = fh->g_spin[(((size_t)p * ns + q) * ns + r) * ns + s];
                    if (fabs(g) < 1e-14) continue;
                    PauliSum AA = pauli_mul(&A[s], &A[r], 1e-15);
                    PauliSum prod = pauli_mul(&CC, &AA, 1e-15);
                    pauli_sum_append(&fh->ham, &prod, 0.5 * g);
                    pauli_sum_free(&prod);
                    pauli_sum_free(&AA);
                }
            pauli_sum_free(&CC);
        }

    for (int j = 0; j < ns; j++) {
        pauli_sum_free(&A[j]);
        pauli_sum_free(&C[j]);
    }
    free(A);
    free(C);

    pauli_sum_merge(&fh->ham, 1e-12);
    return 0;
}

void fermion_free(FermionHam *fh)
{
    pauli_sum_free(&fh->ham);
    free(fh->h_spin);
    free(fh->g_spin);
    memset(fh, 0, sizeof *fh);
}

PauliSum jw_number_operator(int nspin)
{
    /* N = sum_j (1 - Z_j)/2 */
    PauliSum N;
    pauli_sum_init(&N, nspin);
    for (int j = 0; j < nspin; j++) {
        pauli_sum_add(&N, 0, 0, 0.5);
        pauli_sum_add(&N, 0, 1u << j, -0.5);
    }
    pauli_sum_merge(&N, 1e-15);
    return N;
}

/* ---- validation ---- */

static double complex inner(const StateVec *u, const StateVec *v)
{
    double complex s = 0.0;
    for (size_t i = 0; i < u->dim; i++) s += conj(u->a[i]) * v->a[i];
    return s;
}

int fermion_validate(const FermionHam *fh, const ActiveSpace *as, const RHFResult *rhf,
                     double *comm_norm, double *herm_err, double *hf_energy_err, FILE *log)
{
    (void)rhf;
    int ns = fh->nspin;
    int ok = 1;
    int norb = as->nact;

    /* 1. integral symmetries */
    double sym = 0.0;
    for (int s = 0; s < ns; s++)
        for (int t = 0; t < ns; t++) {
            double d = fabs(fh->h_spin[(size_t)s * ns + t] - fh->h_spin[(size_t)t * ns + s]);
            if (d > sym) sym = d;
        }
    for (int p = 0; p < ns; p++)
        for (int q = 0; q < ns; q++)
            for (int r = 0; r < ns; r++)
                for (int s = 0; s < ns; s++) {
                    double v = fh->g_spin[(((size_t)p * ns + q) * ns + r) * ns + s];
                    double w = fh->g_spin[(((size_t)r * ns + s) * ns + p) * ns + q];
                    double d = fabs(v - w);
                    if (d > sym) sym = d;
                }
    /* 2. Hermiticity of the Pauli operator */
    double imax = 0.0;
    int herm = pauli_is_hermitian(&fh->ham, 1e-10, &imax);

    /* 3. number conservation: || [H,N] psi || on a random state */
    PauliSum N = jw_number_operator(ns);
    StateVec psi, a, b, c, d2;
    sv_init(&psi, ns); sv_init(&a, ns); sv_init(&b, ns); sv_init(&c, ns); sv_init(&d2, ns);
    unsigned seed = 12345u;
    double nrm = 0.0;
    for (size_t i = 0; i < psi.dim; i++) {
        seed = seed * 1664525u + 1013904223u;
        double re = ((double)(seed >> 8) / 16777216.0) - 0.5;
        seed = seed * 1664525u + 1013904223u;
        double im = ((double)(seed >> 8) / 16777216.0) - 0.5;
        psi.a[i] = re + I * im;
        nrm += re * re + im * im;
    }
    for (size_t i = 0; i < psi.dim; i++) psi.a[i] /= sqrt(nrm);
    sv_apply_sum(&psi, &N, &a);        /* N psi */
    sv_apply_sum(&a, &fh->ham, &b);    /* H N psi */
    sv_apply_sum(&psi, &fh->ham, &c);  /* H psi */
    sv_apply_sum(&c, &N, &d2);         /* N H psi */
    double cn = 0.0;
    for (size_t i = 0; i < psi.dim; i++) {
        double complex z = b.a[i] - d2.a[i];
        cn += creal(z * conj(z));
    }
    cn = sqrt(cn);

    /* 4. Hermiticity via random vectors: <u|H v> = <H u|v> */
    StateVec u, v;
    sv_init(&u, ns); sv_init(&v, ns);
    for (size_t i = 0; i < u.dim; i++) {
        seed = seed * 1664525u + 1013904223u;
        u.a[i] = ((double)(seed >> 8) / 16777216.0) - 0.5;
        seed = seed * 1664525u + 1013904223u;
        v.a[i] = ((double)(seed >> 8) / 16777216.0) - 0.5 + 0.0 * I;
    }
    sv_apply_sum(&v, &fh->ham, &a);
    sv_apply_sum(&u, &fh->ham, &b);
    double complex lhs = inner(&u, &a);    /* <u|H v> */
    double complex rhs = inner(&b, &v);    /* <H u|v> */
    double hd = cabs(lhs - rhs);

    /* 5. HF determinant: occupied spin orbitals 0..2*nocc-1 */
    int nocc = as->nelec_act / 2;
    uint32_t hf = 0;
    for (int p = 0; p < nocc; p++) hf |= (1u << (2 * p)) | (1u << (2 * p + 1));
    sv_set_basis(&psi, hf);
    sv_apply_sum(&psi, &fh->ham, &a);
    double eh = creal(inner(&psi, &a));
    double e_ref = active_hf_energy(as);
    double hf_err = fabs(eh - e_ref);
    sv_apply_sum(&psi, &N, &a);
    double nh = creal(inner(&psi, &a));

    sv_free(&psi); sv_free(&a); sv_free(&b); sv_free(&c); sv_free(&d2); sv_free(&u); sv_free(&v);
    pauli_sum_free(&N);

    if (comm_norm) *comm_norm = cn;
    if (herm_err) *herm_err = hd;
    if (hf_energy_err) *hf_energy_err = hf_err;

    if (log) {
        fprintf(log, "  integral symmetry violation      : %.3e\n", sym);
        fprintf(log, "  Pauli coefficients max |Im|       : %.3e  (%s)\n", imax, herm ? "Hermitian" : "NOT Hermitian");
        fprintf(log, "  ||[H,N] psi|| (random psi)        : %.3e\n", cn);
        fprintf(log, "  |<u|Hv> - <Hu|v>| (random u,v)     : %.3e\n", hd);
        fprintf(log, "  <HF|N|HF>                         : %.12f  (expected %d)\n", nh, as->nelec_act);
        fprintf(log, "  <HF|H|HF> (Pauli)                 : %.12f\n", eh);
        fprintf(log, "  active-space HF energy (integrals): %.12f\n", e_ref);
        fprintf(log, "  |difference|                      : %.3e\n", hf_err);
        fprintf(log, "  Pauli terms                       : %zu\n", fh->ham.n);
    }
    if (!herm || sym > 1e-10 || cn > 1e-9 || hd > 1e-9 || hf_err > 1e-9 || fabs(nh - as->nelec_act) > 1e-9)
        ok = 0;
    (void)norb;
    return ok ? 0 : -1;
}
