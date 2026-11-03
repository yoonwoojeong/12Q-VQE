/*
 * Restricted Hartree-Fock in the AO basis.
 * Löwdin orthogonalisation, LAPACK dsyev, BLAS dgemm, Pulay DIIS.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DIIS_MAX 8

/* C = A*B (or with transposes), all n x n row-major */
static void mm(int n, const double *A, int ta, const double *B, int tb, double *C)
{
    cblas_dgemm(CblasRowMajor, ta ? CblasTrans : CblasNoTrans, tb ? CblasTrans : CblasNoTrans,
                n, n, n, 1.0, A, n, B, n, 0.0, C, n);
}

static void build_fock(int n, const double *H, const double *P, const double *eri, double *F)
{
    /* G_pq = sum_rs P_rs [ (pq|rs) - 1/2 (pr|qs) ] */
    for (int p = 0; p < n; p++)
        for (int q = 0; q < n; q++) {
            double g = 0.0;
            for (int r = 0; r < n; r++) {
                const double *erir = eri + (((size_t)p * n + q) * n + r) * n; /* (pq|r.) */
                for (int s = 0; s < n; s++) {
                    double Prs = P[(size_t)r * n + s];
                    if (Prs == 0.0) continue;
                    double jpq = erir[s];
                    double kprqs = eri[(((size_t)p * n + r) * n + q) * n + s];
                    g += Prs * (jpq - 0.5 * kprqs);
                }
            }
            F[(size_t)p * n + q] = H[(size_t)p * n + q] + g;
        }
}

/* Solve F C = eps S C via the orthogonal basis X; returns C (AO x MO) */
static void diag_fock(int n, const double *F, const double *X, double *Cout, double *eps)
{
    double *Fp = qc_calloc((size_t)n * n, sizeof(double));
    double *tmp = qc_calloc((size_t)n * n, sizeof(double));
    mm(n, X, 1, F, 0, tmp);      /* X^T F */
    mm(n, tmp, 0, X, 0, Fp);     /* X^T F X */
    /* symmetrise against round-off */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++) {
            double v = 0.5 * (Fp[(size_t)i * n + j] + Fp[(size_t)j * n + i]);
            Fp[(size_t)i * n + j] = Fp[(size_t)j * n + i] = v;
        }
    double *w = qc_calloc((size_t)n, sizeof(double));
    if (LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, Fp, n, w) != 0) qc_die("dsyev failed in RHF");
    memcpy(eps, w, sizeof(double) * (size_t)n);
    mm(n, X, 0, Fp, 0, Cout);   /* C = X C' */
    free(w);
    free(Fp);
    free(tmp);
}

static void density_from_C(int n, const double *C, int nocc, double *P)
{
    memset(P, 0, sizeof(double) * (size_t)n * n);
    for (int m = 0; m < n; m++)
        for (int v = 0; v < n; v++) {
            double s = 0.0;
            for (int i = 0; i < nocc; i++) s += C[(size_t)m * n + i] * C[(size_t)v * n + i];
            P[(size_t)m * n + v] = 2.0 * s;
        }
}

int rhf_run(RHFResult *r, const AOInts *ints, const Molecule *mol, int maxiter, double econv, double dconv, int verbose)
{
    int n = ints->nao;
    memset(r, 0, sizeof *r);
    r->nao = n;
    r->enuc = mol->enuc;
    if (mol->nelec % 2 != 0 || mol->nelec <= 0) qc_die("RHF requires a closed-shell electron count (nelec=%d)", mol->nelec);
    int nocc = mol->nelec / 2;
    if (nocc > n) qc_die("RHF: more occupied orbitals than basis functions");
    r->nocc = nocc;

    /* X = S^{-1/2} */
    double *Us = qc_calloc((size_t)n * n, sizeof(double));
    double *sv = qc_calloc((size_t)n, sizeof(double));
    memcpy(Us, ints->S, sizeof(double) * (size_t)n * n);
    if (LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, Us, n, sv) != 0) qc_die("dsyev failed on overlap");
    double *X = qc_calloc((size_t)n * n, sizeof(double));
    for (int i = 0; i < n; i++) {
        if (sv[i] < 1e-10) qc_die("RHF: overlap matrix is near-singular (eigenvalue %.3e)", sv[i]);
        double f = 1.0 / sqrt(sv[i]);
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
                X[(size_t)a * n + b] += Us[(size_t)a * n + i] * f * Us[(size_t)b * n + i];
    }
    free(sv);
    free(Us);

    r->C = qc_calloc((size_t)n * n, sizeof(double));
    r->eps = qc_calloc((size_t)n, sizeof(double));
    r->P = qc_calloc((size_t)n * n, sizeof(double));
    double *F = qc_calloc((size_t)n * n, sizeof(double));
    double *Pold = qc_calloc((size_t)n * n, sizeof(double));
    double *FPS = qc_calloc((size_t)n * n, sizeof(double));
    double *tmp1 = qc_calloc((size_t)n * n, sizeof(double));
    double *tmp2 = qc_calloc((size_t)n * n, sizeof(double));
    double *err = qc_calloc((size_t)n * n, sizeof(double));

    /* DIIS storage */
    double *Fh[DIIS_MAX], *Eh[DIIS_MAX];
    int nd = 0;
    for (int k = 0; k < DIIS_MAX; k++) { Fh[k] = NULL; Eh[k] = NULL; }

    /* core-Hamiltonian guess */
    diag_fock(n, ints->H, X, r->C, r->eps);
    density_from_C(n, r->C, nocc, r->P);

    double Eold = 0.0;
    int conv = 0;
    double Etot = 0.0;
    double dP_rms = 0.0;
    for (int it = 1; it <= maxiter; it++) {
        build_fock(n, ints->H, r->P, ints->eri, F);

        /* energy E = 1/2 sum P (H+F) + Enuc */
        double E = 0.0;
        for (size_t k = 0; k < (size_t)n * n; k++) E += 0.5 * r->P[k] * (ints->H[k] + F[k]);
        Etot = E + mol->enuc;

        /* DIIS error e = X^T (F P S - S P F) X */
        mm(n, F, 0, r->P, 0, tmp1);
        mm(n, tmp1, 0, ints->S, 0, FPS);            /* F P S */
        mm(n, ints->S, 0, r->P, 0, tmp1);
        mm(n, tmp1, 0, F, 0, tmp2);                 /* S P F */
        for (size_t k = 0; k < (size_t)n * n; k++) tmp1[k] = FPS[k] - tmp2[k];
        mm(n, X, 1, tmp1, 0, tmp2);
        mm(n, tmp2, 0, X, 0, err);
        double emax = 0.0;
        for (size_t k = 0; k < (size_t)n * n; k++) if (fabs(err[k]) > emax) emax = fabs(err[k]);

        /* DIIS extrapolation */
        if (nd == DIIS_MAX) {
            free(Fh[0]); free(Eh[0]);
            memmove(Fh, Fh + 1, sizeof(double *) * (DIIS_MAX - 1));
            memmove(Eh, Eh + 1, sizeof(double *) * (DIIS_MAX - 1));
            nd--;
        }
        Fh[nd] = qc_calloc((size_t)n * n, sizeof(double));
        Eh[nd] = qc_calloc((size_t)n * n, sizeof(double));
        memcpy(Fh[nd], F, sizeof(double) * (size_t)n * n);
        memcpy(Eh[nd], err, sizeof(double) * (size_t)n * n);
        nd++;

        double *Fuse = F;
        double *Fdiis = NULL;
        if (nd >= 2) {
            int m = nd + 1;
            double *B = qc_calloc((size_t)m * m, sizeof(double));
            double *rhs = qc_calloc((size_t)m, sizeof(double));
            int *ipiv = qc_calloc((size_t)m, sizeof(int));
            for (int a = 0; a < nd; a++)
                for (int b = 0; b < nd; b++) {
                    double s = 0.0;
                    for (size_t k = 0; k < (size_t)n * n; k++) s += Eh[a][k] * Eh[b][k];
                    B[(size_t)a * m + b] = s;
                }
            for (int a = 0; a < nd; a++) {
                B[(size_t)a * m + nd] = -1.0;
                B[(size_t)nd * m + a] = -1.0;
            }
            B[(size_t)nd * m + nd] = 0.0;
            rhs[nd] = -1.0;
            if (LAPACKE_dgesv(LAPACK_ROW_MAJOR, m, 1, B, m, ipiv, rhs, 1) == 0) {
                Fdiis = qc_calloc((size_t)n * n, sizeof(double));
                for (int a = 0; a < nd; a++)
                    for (size_t k = 0; k < (size_t)n * n; k++) Fdiis[k] += rhs[a] * Fh[a][k];
                Fuse = Fdiis;
            }
            free(B); free(rhs); free(ipiv);
        }

        /* diagonalise, build new density */
        diag_fock(n, Fuse, X, r->C, r->eps);
        memcpy(Pold, r->P, sizeof(double) * (size_t)n * n);
        density_from_C(n, r->C, nocc, r->P);
        double s2 = 0.0;
        for (size_t k = 0; k < (size_t)n * n; k++) {
            double d = r->P[k] - Pold[k];
            s2 += d * d;
        }
        dP_rms = sqrt(s2 / ((double)n * n));
        double dE = Etot - Eold;
        if (verbose)
            printf("  RHF iter %3d  E = %.12f  dE = %+.3e  |[F,D]|max = %.3e  dP_rms = %.3e\n", it, Etot, dE, emax, dP_rms);
        r->iterations = it;
        r->delta_e = dE;
        r->rms_density = dP_rms;
        Eold = Etot;
        free(Fdiis);
        if (it > 1 && fabs(dE) < econv && dP_rms < dconv) { conv = 1; break; }
    }

    /* final energy with the converged density */
    build_fock(n, ints->H, r->P, ints->eri, F);
    double E = 0.0;
    for (size_t k = 0; k < (size_t)n * n; k++) E += 0.5 * r->P[k] * (ints->H[k] + F[k]);
    r->energy = E + mol->enuc;

    for (int k = 0; k < DIIS_MAX; k++) { free(Fh[k]); free(Eh[k]); }
    free(F); free(Pold); free(FPS); free(tmp1); free(tmp2); free(err); free(X);
    if (!conv) qc_warn("RHF did not converge in %d iterations", maxiter);
    (void)Etot;
    return conv ? 0 : -1;
}

void rhf_free(RHFResult *r)
{
    free(r->C);
    free(r->eps);
    free(r->P);
    memset(r, 0, sizeof *r);
}
