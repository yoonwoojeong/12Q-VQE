/*
 * Active-space construction: AO -> MO integral transformation for a window of
 * spatial orbitals [first, first+norb). Doubly occupied orbitals below the
 * window are inactive and are folded into an effective one-body operator and a
 * constant energy (E_core).
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* out[a][b][c][d] = sum (mu nu|la si) Ca[mu][a] Cb[nu][b] Cc[la][c] Cd[si][d] (chemist) */
static double *ao2mo4(int n, const double *eri, const double *Ca, int na, const double *Cb, int nb,
                      const double *Cc, int nc, const double *Cd, int nd)
{
    size_t n3 = (size_t)n * n * n;
    /* step 1: contract first index */
    double *T1 = qc_calloc((size_t)na * n3, sizeof(double));
    for (int a = 0; a < na; a++)
        for (size_t mu = 0; mu < (size_t)n; mu++) {
            double c = Ca[mu * na + a];
            if (c == 0.0) continue;
            const double *src = eri + mu * n3;
            double *dst = T1 + (size_t)a * n3;
            for (size_t k = 0; k < n3; k++) dst[k] += c * src[k];
        }
    /* step 2: second index nu -> b: T2[a][b][la][si] */
    double *T2 = qc_calloc((size_t)na * nb * n * n, sizeof(double));
    for (int a = 0; a < na; a++)
        for (int b = 0; b < nb; b++)
            for (int nu = 0; nu < n; nu++) {
                double c = Cb[(size_t)nu * nb + b];
                if (c == 0.0) continue;
                const double *src = T1 + (size_t)a * n3 + (size_t)nu * n * n;
                double *dst = T2 + ((size_t)a * nb + b) * n * n;
                for (size_t k = 0; k < (size_t)n * n; k++) dst[k] += c * src[k];
            }
    free(T1);
    /* step 3: la -> c : T3[a][b][c][si] */
    double *T3 = qc_calloc((size_t)na * nb * nc * n, sizeof(double));
    for (int a = 0; a < na; a++)
        for (int b = 0; b < nb; b++)
            for (int c = 0; c < nc; c++)
                for (int la = 0; la < n; la++) {
                    double cc = Cc[(size_t)la * nc + c];
                    if (cc == 0.0) continue;
                    const double *src = T2 + (((size_t)a * nb + b) * n + la) * n;
                    double *dst = T3 + (((size_t)a * nb + b) * nc + c) * n;
                    for (int si = 0; si < n; si++) dst[si] += cc * src[si];
                }
    free(T2);
    /* step 4: si -> d : out[a][b][c][d] */
    double *out = qc_calloc((size_t)na * nb * nc * nd, sizeof(double));
    for (int a = 0; a < na; a++)
        for (int b = 0; b < nb; b++)
            for (int c = 0; c < nc; c++)
                for (int d = 0; d < nd; d++) {
                    double s = 0.0;
                    const double *src = T3 + (((size_t)a * nb + b) * nc + c) * n;
                    for (int si = 0; si < n; si++) s += src[si] * Cd[(size_t)si * nd + d];
                    out[(((size_t)a * nb + b) * nc + c) * nd + d] = s;
                }
    free(T3);
    return out;
}

/* one-electron MO matrix: out[p][q] = sum Ca[mu][p] Cb[nu][q] M[mu][nu] */
static void mo1e(int n, const double *M, const double *Ca, int na, const double *Cb, int nb, double *out)
{
    double *tmp = qc_calloc((size_t)n * nb, sizeof(double));
    for (int mu = 0; mu < n; mu++)
        for (int b = 0; b < nb; b++) {
            double s = 0.0;
            for (int nu = 0; nu < n; nu++) s += M[(size_t)mu * n + nu] * Cb[(size_t)nu * nb + b];
            tmp[(size_t)mu * nb + b] = s;
        }
    for (int a = 0; a < na; a++)
        for (int b = 0; b < nb; b++) {
            double s = 0.0;
            for (int mu = 0; mu < n; mu++) s += Ca[(size_t)mu * na + a] * tmp[(size_t)mu * nb + b];
            out[(size_t)a * nb + b] = s;
        }
    free(tmp);
}

int active_space_build(ActiveSpace *as, const AOInts *ints, const RHFResult *rhf,
                       const Molecule *mol, int nelec_act, int norb, int first)
{
    int n = ints->nao;
    memset(as, 0, sizeof *as);
    as->nao = n;
    as->nact = norb;
    as->first = first;
    as->nelec_act = nelec_act;
    if (first < 0 || norb <= 0 || first + norb > n) qc_die("active space window [%d,%d) is outside the %d MOs", first, first + norb, n);
    if (mol->nelec - nelec_act != 2 * first)
        qc_die("active space: %d total electrons, %d active, so the window must start at MO %d (inactive doubly occupied MOs below it)",
               mol->nelec, nelec_act, (mol->nelec - nelec_act) / 2);
    if (nelec_act < 0 || nelec_act > 2 * norb) qc_die("active electron count %d does not fit in %d orbitals", nelec_act, norb);
    if (nelec_act % 2 != 0) qc_die("only closed-shell active spaces (even electron count) are supported");
    int ncore = first;
    as->ncore = ncore;

    /* MO coefficient slices from the RHF solution (columns are MOs) */
    double *Cact = qc_calloc((size_t)n * norb, sizeof(double));
    double *Ccore = qc_calloc((size_t)n * (ncore ? ncore : 1), sizeof(double));
    for (int mu = 0; mu < n; mu++) {
        for (int k = 0; k < norb; k++) Cact[(size_t)mu * norb + k] = rhf->C[(size_t)mu * n + first + k];
        for (int i = 0; i < ncore; i++) Ccore[(size_t)mu * ncore + i] = rhf->C[(size_t)mu * n + i];
    }
    as->C_act = Cact;
    as->eps_act = qc_calloc((size_t)norb, sizeof(double));
    for (int k = 0; k < norb; k++) as->eps_act[k] = rhf->eps[first + k];

    /* one-electron MO integrals */
    double *hmo_act = qc_calloc((size_t)norb * norb, sizeof(double));
    mo1e(n, ints->H, Cact, norb, Cact, norb, hmo_act);
    double *hmo_core = qc_calloc((size_t)(ncore ? ncore : 1) * (ncore ? ncore : 1), sizeof(double));
    if (ncore) mo1e(n, ints->H, Ccore, ncore, Ccore, ncore, hmo_core);

    /* core energy and effective one-body operator */
    double ecore = mol->enuc;
    for (int i = 0; i < ncore; i++) ecore += 2.0 * hmo_core[(size_t)i * ncore + i];
    as->h = qc_calloc((size_t)norb * norb, sizeof(double));
    memcpy(as->h, hmo_act, sizeof(double) * (size_t)norb * norb);
    if (ncore) {
        double *cc = ao2mo4(n, ints->eri, Ccore, ncore, Ccore, ncore, Ccore, ncore, Ccore, ncore);   /* (ii|jj) */
        for (int i = 0; i < ncore; i++)
            for (int j = 0; j < ncore; j++) {
                double J = cc[(((size_t)i * ncore + i) * ncore + j) * ncore + j];
                double K = cc[(((size_t)i * ncore + j) * ncore + j) * ncore + i];   /* (ij|ji) */
                ecore += 2.0 * J - K;
            }
        free(cc);
        /* effective one-body: h_pq += sum_c [2 (pq|cc) - (pc|cq)] */
        double *pqcc = ao2mo4(n, ints->eri, Cact, norb, Cact, norb, Ccore, ncore, Ccore, ncore);
        double *pccq = ao2mo4(n, ints->eri, Cact, norb, Ccore, ncore, Ccore, ncore, Cact, norb);
        for (int p = 0; p < norb; p++)
            for (int q = 0; q < norb; q++)
                for (int c = 0; c < ncore; c++) {
                    double J = pqcc[(((size_t)p * norb + q) * ncore + c) * ncore + c];
                    double K = pccq[(((size_t)p * ncore + c) * ncore + c) * norb + q];
                    as->h[(size_t)p * norb + q] += 2.0 * J - K;
                }
        free(pqcc);
        free(pccq);
    }
    as->ecore = ecore;
    free(hmo_act);
    free(hmo_core);

    /* active two-electron integrals (pq|rs) chemist order */
    as->g = ao2mo4(n, ints->eri, Cact, norb, Cact, norb, Cact, norb, Cact, norb);
    free(Ccore);
    return 0;
}

double active_hf_energy(const ActiveSpace *as)
{
    int na = as->nact;
    int nocc = as->nelec_act / 2;
    double e = as->ecore;
    for (int i = 0; i < nocc; i++) e += 2.0 * as->h[(size_t)i * na + i];
    for (int i = 0; i < nocc; i++)
        for (int j = 0; j < nocc; j++) {
            double J = as->g[(((size_t)i * na + i) * na + j) * na + j];
            double K = as->g[(((size_t)i * na + j) * na + j) * na + i];
            e += 2.0 * J - K;
        }
    return e;
}

void active_space_free(ActiveSpace *as)
{
    free(as->h);
    free(as->g);
    free(as->C_act);
    free(as->eps_act);
    memset(as, 0, sizeof *as);
}
