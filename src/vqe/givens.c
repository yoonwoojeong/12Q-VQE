/*
 * Shallow state preparation for a two-electron (one alpha, one beta) state.
 *
 * Any such state can be written as
 *     |Psi> = sum_pq C_pq a+_{p alpha} a+_{q beta} |0>,
 * and the SVD  C = U S V^T  gives its natural-orbital (paired) form
 *     |Psi> = sum_k s_k b+_{k alpha} c+_{k beta} |0>,
 *     b+_k = sum_p U_pk a+_{p alpha},  c+_k = sum_q V_qk a+_{q beta}.
 * The circuit
 *   1. puts the pair amplitudes s_k on modes k = 0..r-1 with a chain of pair rotations
 *      exp(phi (P+_{k+1} P_k - h.c.)),  P+_k = a+_{k alpha} a+_{k beta},
 *   2. rotates alpha and beta orbitals with adjacent Givens rotations
 *      exp(theta (a+_i a_j - a+_j a_i)) so that mode k becomes natural orbital k.
 * Each pair rotation acts on 4 consecutive qubits and each Givens rotation on 3
 * (spin-interleaved Jordan-Wigner ordering), so on a path layout no SWAPs are needed.
 * The circuit is checked against the target on the statevector simulator.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <lapacke.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NMAX (QC_MAX_QUBIT_SV / 2)

typedef struct {
    int i, j;
    double theta;
} GivensOp;

/* Zero the first r columns of the orthogonal matrix M (n x n, row-major) below the diagonal
   with adjacent-row rotations R(i,j,theta) = [[c, s], [-s, c]]. Records the rotations and
   returns the signs d_k = M[k][k] that remain on the diagonal. */
static int eliminate(double *M, int n, int r, GivensOp *ops, double *d)
{
    int nop = 0;
    for (int c = 0; c < r; c++) {
        for (int row = n - 1; row > c; row--) {
            int i = row - 1, j = row;
            double a = M[i * n + c], b = M[j * n + c];
            if (fabs(b) < 1e-14) continue;
            double th = atan2(b, a);
            double cs = cos(th), sn = sin(th);
            for (int col = 0; col < n; col++) {
                double mi = M[i * n + col], mj = M[j * n + col];
                M[i * n + col] = cs * mi + sn * mj;
                M[j * n + col] = -sn * mi + cs * mj;
            }
            ops[nop].i = i;
            ops[nop].j = j;
            ops[nop].theta = th;
            nop++;
        }
        d[c] = M[c * n + c] >= 0 ? 1.0 : -1.0;
    }
    return nop;
}

static void append_generator(Circuit *c, int nc, const int *cre, int na, const int *ann, double angle, int nq)
{
    GenTerm *t = NULL;
    int nt = 0;
    if (fermion_generator(nc, cre, na, ann, &t, &nt) != 0) qc_die("givens: generator is not anti-Hermitian");
    for (int k = 0; k < nt; k++) circ_append_pauli_rotation(c, t[k].x, t[k].z, angle * t[k].d, nq);
    free(t);
}

int givens_state_circuit(const StateVec *target, int nspin, double trunc, Circuit *c, GivensInfo *info)
{
    int n = nspin / 2;
    if (nspin % 2 || n > NMAX || target->nq != nspin) qc_die("givens: bad register size");
    memset(info, 0, sizeof *info);

    /* coefficient matrix C_pq of a+_{p alpha} a+_{q beta}|0>; JW basis states are ascending
       products, so the amplitude of bits {2p, 2q+1} carries a sign when 2q+1 < 2p */
    double C[NMAX * NMAX];
    double inside = 0.0, total = 0.0, maxim = 0.0;
    for (size_t idx = 0; idx < target->dim; idx++) total += creal(target->a[idx] * conj(target->a[idx]));
    for (int p = 0; p < n; p++)
        for (int q = 0; q < n; q++) {
            size_t idx = ((size_t)1 << (2 * p)) | ((size_t)1 << (2 * q + 1));
            double complex a = target->a[idx];
            if (fabs(cimag(a)) > maxim) maxim = fabs(cimag(a));
            C[p * n + q] = (p <= q ? 1.0 : -1.0) * creal(a);
            inside += creal(a * conj(a));
        }
    info->leakage = total > 0 ? 1.0 - inside / total : 1.0;
    if (info->leakage > 1e-8) qc_die("givens: target has weight %.3e outside the one-alpha/one-beta sector", info->leakage);
    if (maxim > 1e-8) qc_die("givens: target amplitudes are not real (max |Im| %.3e)", maxim);

    /* SVD  C = U S V^T  (singular values descending) */
    double A[NMAX * NMAX], S[NMAX], U[NMAX * NMAX], VT[NMAX * NMAX], superb[NMAX];
    memcpy(A, C, sizeof(double) * (size_t)(n * n));
    if (LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', n, n, A, n, S, U, n, VT, n, superb) != 0)
        qc_die("givens: SVD failed");
    double norm2 = 0.0;
    for (int k = 0; k < n; k++) norm2 += S[k] * S[k];
    int r = 0;
    double kept = 0.0;
    while (r < n && kept < (1.0 - trunc) * norm2) { kept += S[r] * S[r]; r++; }
    if (r == 0) r = 1;
    info->r = r;
    info->kept_weight = kept / norm2;
    for (int k = 0; k < n; k++) info->sv[k] = S[k] / sqrt(norm2);

    /* Givens decompositions of U (alpha) and V (beta) for the first r natural orbitals */
    double V[NMAX * NMAX];
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++) V[a * n + b] = VT[b * n + a];
    GivensOp opa[NMAX * NMAX], opb[NMAX * NMAX];
    double da[NMAX], db[NMAX];
    int na = eliminate(U, n, r, opa, da);
    int nb = eliminate(V, n, r, opb, db);

    /* signed, normalised pair amplitudes t_k = s_k d^alpha_k d^beta_k */
    double t[NMAX], tn = 0.0;
    for (int k = 0; k < r; k++) { t[k] = S[k] * da[k] * db[k]; tn += t[k] * t[k]; }
    for (int k = 0; k < r; k++) t[k] /= sqrt(tn);

    circ_init(c, nspin);
    circ_add(c, G_X, 0, -1, 0);   /* |0 alpha, 0 beta> */
    circ_add(c, G_X, 1, -1, 0);
    if (t[0] < 0 && r == 1) { /* global sign only */ }

    /* 1. pair rotations k -> k+1 with amplitudes a_k = t_k */
    for (int k = 0; k + 1 < r; k++) {
        double phi;
        if (k < r - 2) {
            double rest = 0.0;
            for (int j = k + 1; j < r; j++) rest += t[j] * t[j];
            phi = atan2(sqrt(rest), t[k]);
        } else {
            phi = atan2(t[k + 1], t[k]);
        }
        int cre[2] = {2 * (k + 1), 2 * (k + 1) + 1};
        int ann[2] = {2 * k + 1, 2 * k};
        append_generator(c, 2, cre, 2, ann, phi, nspin);
        info->n_pair_rot++;
    }
    /* overall sign of the chain when t_0 < 0 for r > 1 is handled by atan2 above, except the
       first amplitude: cos(phi_0) carries the sign of t_0 */

    /* 2. orbital rotations: U = G_1^T ... G_L^T D, applied as G_L^T first, each R(i,j,-theta) */
    for (int m = na - 1; m >= 0; m--) {
        int cre[1] = {2 * opa[m].i}, ann[1] = {2 * opa[m].j};
        append_generator(c, 1, cre, 1, ann, -opa[m].theta, nspin);
        info->n_givens++;
    }
    for (int m = nb - 1; m >= 0; m--) {
        int cre[1] = {2 * opb[m].i + 1}, ann[1] = {2 * opb[m].j + 1};
        append_generator(c, 1, cre, 1, ann, -opb[m].theta, nspin);
        info->n_givens++;
    }

    /* verification on the simulator */
    StateVec sv;
    sv_init(&sv, nspin);
    circ_simulate(c, &sv);
    double complex ov = 0.0;
    for (size_t idx = 0; idx < sv.dim; idx++) ov += conj(target->a[idx]) * sv.a[idx];
    info->fidelity = creal(ov * conj(ov)) / total;
    sv_free(&sv);
    return 0;
}
