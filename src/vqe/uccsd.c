/*
 * Spin-conserving UCCSD ansatz from the Hartree-Fock determinant.
 *
 * Excitation operators (spin orbitals, JW-mapped):
 *   singles  T = a+_a a_i             (same spin)
 *   doubles  T = a+_a a+_b a_j a_i    (total Sz conserved)
 * Generator G = T - T+ = sum_s i d_s P_s  (P_s canonical Hermitian Pauli).
 * exp(theta G) is approximated by the first-order Trotter product
 *   prod_s exp(i theta d_s P_s),  exact when the strings of one excitation commute.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* product of ladder operators given as (orbital, dagger) pairs, left to right */
static PauliSum op_product(int n, const int *orb, const int *dag)
{
    PauliSum acc;
    pauli_sum_init(&acc, 0);
    PauliSum one;
    pauli_sum_init(&one, 0);
    pauli_sum_add(&one, 0, 0, 1.0);
    /* start from identity */
    pauli_sum_init(&acc, 0);
    pauli_sum_add(&acc, 0, 0, 1.0);
    for (int k = 0; k < n; k++) {
        PauliSum L;
        fermion_ladder(orb[k], dag[k], &L);
        PauliSum prod = pauli_mul(&acc, &L, 1e-15);
        pauli_sum_free(&acc);
        pauli_sum_free(&L);
        acc = prod;
    }
    pauli_sum_free(&one);
    return acc;
}

/* Pauli terms of G = T - T+ for T = a+_{cre[0]} ... a+_{cre[nc-1]} a_{ann[0]} ... a_{ann[na-1]}.
   G is anti-Hermitian, so G = sum_s i d_s P_s; returns the (x, z, d) list. */
int fermion_generator(int nc, const int *cre, int na, const int *ann, GenTerm **out, int *nout)
{
    int orb[8], dag[8], orbd[8], dagd[8];
    int n = nc + na;
    if (n > 8) return -1;
    for (int k = 0; k < nc; k++) { orb[k] = cre[k]; dag[k] = 1; }
    for (int k = 0; k < na; k++) { orb[nc + k] = ann[k]; dag[nc + k] = 0; }
    /* T+ : reverse order, flip daggers */
    for (int k = 0; k < n; k++) { orbd[k] = orb[n - 1 - k]; dagd[k] = !dag[n - 1 - k]; }
    PauliSum T = op_product(n, orb, dag);
    PauliSum Td = op_product(n, orbd, dagd);
    PauliSum G;
    pauli_sum_init(&G, 0);
    pauli_sum_append(&G, &T, 1.0);
    pauli_sum_append(&G, &Td, -1.0);
    pauli_sum_merge(&G, 1e-14);
    pauli_sum_free(&T);
    pauli_sum_free(&Td);
    *out = qc_calloc(G.n ? G.n : 1, sizeof(GenTerm));
    *nout = 0;
    for (size_t t = 0; t < G.n; t++) {
        if (fabs(creal(G.t[t].c)) > 1e-10) { pauli_sum_free(&G); free(*out); *out = NULL; return -1; }
        (*out)[*nout].x = G.t[t].x;
        (*out)[*nout].z = G.t[t].z;
        (*out)[*nout].d = cimag(G.t[t].c);
        (*nout)++;
    }
    pauli_sum_free(&G);
    return 0;
}

static int spin_of(int s) { return s & 1; }

int uccsd_build(UCCSD *u, int nspin, int nelec)
{
    memset(u, 0, sizeof *u);
    u->nq = nspin;
    u->nelec = nelec;
    if (nelec % 2) qc_die("UCCSD: closed-shell reference required");
    int nocc_sp = nelec;                          /* occupied spin orbitals 0..nelec-1 (2 per spatial, pairs) */
    /* occupied spin orbitals: 2p and 2p+1 for p < nelec/2 */
    int occ[QC_MAX_QUBIT_SV], vir[QC_MAX_QUBIT_SV];
    int no = 0, nv = 0;
    for (int s = 0; s < nspin; s++) {
        int p = s >> 1;
        if (p < nelec / 2) occ[no++] = s;
        else vir[nv++] = s;
    }
    (void)nocc_sp;
    u->hf_mask = 0;
    for (int k = 0; k < no; k++) u->hf_mask |= 1u << occ[k];
    u->hf_index = u->hf_mask;

    int cap = no * nv + no * (no - 1) / 2 * nv * (nv - 1) / 2 + 4;
    u->exc = qc_calloc((size_t)cap, sizeof(Excitation));
    int ne = 0;

    /* singles */
    for (int i = 0; i < no; i++)
        for (int a = 0; a < nv; a++) {
            if (spin_of(occ[i]) != spin_of(vir[a])) continue;
            Excitation *e = &u->exc[ne++];
            e->kind = 1;
            e->occ[0] = occ[i];
            e->vir[0] = vir[a];
            int orb[2] = {vir[a], occ[i]};
            int dag[2] = {1, 0};
            int ord_d[2] = {occ[i], vir[a]};
            int dagd[2] = {1, 0};
            (void)ord_d; (void)dagd;
            PauliSum T = op_product(2, orb, dag);
            int orbd[2] = {occ[i], vir[a]};
            int dagTd[2] = {1, 0};
            PauliSum Td = op_product(2, orbd, dagTd);
            PauliSum G;
            pauli_sum_init(&G, nspin);
            pauli_sum_append(&G, &T, 1.0);
            pauli_sum_append(&G, &Td, -1.0);
            pauli_sum_merge(&G, 1e-14);
            /* Td computed as a+_i a_a: fix orientation below */
            pauli_sum_free(&T);
            pauli_sum_free(&Td);
            /* store */
            e->terms = qc_calloc(G.n ? G.n : 1, sizeof(GenTerm));
            e->nterm = 0;
            for (size_t t = 0; t < G.n; t++) {
                double complex c = G.t[t].c;
                if (fabs(creal(c)) > 1e-10) qc_die("UCCSD: generator not anti-Hermitian (Re c = %g)", creal(c));
                e->terms[e->nterm].x = G.t[t].x;
                e->terms[e->nterm].z = G.t[t].z;
                e->terms[e->nterm].d = cimag(c);
                e->nterm++;
            }
            pauli_sum_free(&G);
        }

    /* doubles: i<j in occ, a<b in vir, same spin multiset */
    for (int i = 0; i < no; i++)
        for (int j = i + 1; j < no; j++)
            for (int a = 0; a < nv; a++)
                for (int b = a + 1; b < nv; b++) {
                    int si = spin_of(occ[i]), sj = spin_of(occ[j]);
                    int sa = spin_of(vir[a]), sb = spin_of(vir[b]);
                    if (si + sj != sa + sb) continue;   /* Sz conservation: equal beta count */
                    Excitation *e = &u->exc[ne++];
                    e->kind = 2;
                    e->occ[0] = occ[i]; e->occ[1] = occ[j];
                    e->vir[0] = vir[a]; e->vir[1] = vir[b];
                    /* T = a+_a a+_b a_j a_i */
                    int orb[4] = {vir[a], vir[b], occ[j], occ[i]};
                    int dag[4] = {1, 1, 0, 0};
                    PauliSum T = op_product(4, orb, dag);
                    /* T+ = a+_i a+_j a_b a_a */
                    int orbd[4] = {occ[i], occ[j], vir[b], vir[a]};
                    int dagd[4] = {1, 1, 0, 0};
                    PauliSum Td = op_product(4, orbd, dagd);
                    PauliSum G;
                    pauli_sum_init(&G, nspin);
                    pauli_sum_append(&G, &T, 1.0);
                    pauli_sum_append(&G, &Td, -1.0);
                    pauli_sum_merge(&G, 1e-14);
                    pauli_sum_free(&T);
                    pauli_sum_free(&Td);
                    e->terms = qc_calloc(G.n ? G.n : 1, sizeof(GenTerm));
                    e->nterm = 0;
                    for (size_t t = 0; t < G.n; t++) {
                        double complex c = G.t[t].c;
                        if (fabs(creal(c)) > 1e-10) qc_die("UCCSD: generator not anti-Hermitian");
                        e->terms[e->nterm].x = G.t[t].x;
                        e->terms[e->nterm].z = G.t[t].z;
                        e->terms[e->nterm].d = cimag(c);
                        e->nterm++;
                    }
                    pauli_sum_free(&G);
                }
    u->nexc = ne;
    return 0;
}

void uccsd_free(UCCSD *u)
{
    for (int k = 0; k < u->nexc; k++) free(u->exc[k].terms);
    free(u->exc);
    memset(u, 0, sizeof *u);
}

void uccsd_prepare(const UCCSD *u, const double *theta, StateVec *sv)
{
    sv_set_basis(sv, u->hf_index);
    for (int k = 0; k < u->nexc; k++) {
        const Excitation *e = &u->exc[k];
        for (int t = 0; t < e->nterm; t++)
            sv_apply_pauli_rotation(sv, e->terms[t].x, e->terms[t].z, theta[k] * e->terms[t].d);
    }
}

void uccsd_circuit(const UCCSD *u, const double *theta, Circuit *c)
{
    circ_init(c, u->nq);
    for (int q = 0; q < u->nq; q++)
        if ((u->hf_mask >> q) & 1) circ_add(c, G_X, q, -1, 0);
    for (int k = 0; k < u->nexc; k++) {
        const Excitation *e = &u->exc[k];
        for (int t = 0; t < e->nterm; t++)
            circ_append_pauli_rotation(c, e->terms[t].x, e->terms[t].z, theta[k] * e->terms[t].d, u->nq);
    }
}
