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
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        if (dets[mid] == m) return mid;
        if (dets[mid] < m) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;
}

static int cmp_u32(const void *a, const void *b)
{
    uint32_t x = *(const uint32_t *)a;
    uint32_t y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

int fci_solve_subspace(FCIResult *r, const FermionHam *fh, const ActiveSpace *as,
                       const uint32_t *dets, int ndet)
{
    memset(r, 0, sizeof *r);
    if (ndet <= 0) return -1;
    int ns = fh->nspin;

    /* Copy, sort, and deduplicate determinants */
    uint32_t *uniq = qc_calloc((size_t)ndet, sizeof(uint32_t));
    memcpy(uniq, dets, (size_t)ndet * sizeof(uint32_t));
    qsort(uniq, (size_t)ndet, sizeof(uint32_t), cmp_u32);
    int nd = 0;
    for (int i = 0; i < ndet; i++) {
        if (i == 0 || uniq[i] != uniq[i - 1]) uniq[nd++] = uniq[i];
    }
    r->dets = uniq;
    r->ndet = nd;

    double *H = qc_calloc((size_t)nd * nd, sizeof(double));

    #pragma omp parallel for
    for (int J = 0; J < nd; J++) {
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

    double *w = qc_calloc((size_t)nd, sizeof(double));
    if (LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nd, H, nd, w) != 0) qc_die("dsyev failed in subspace FCI");
    r->eigvals = w;
    r->ground = w[0];
    r->ground_vec = qc_calloc((size_t)nd, sizeof(double));
    for (int i = 0; i < nd; i++) r->ground_vec[i] = H[(size_t)i * nd + 0];
    free(H);
    return 0;
}

int fci_solve(FCIResult *r, const FermionHam *fh, const ActiveSpace *as, int nelec)
{
    int ns = fh->nspin;
    int nd = 0;
    uint32_t *all = NULL;
    enumerate_dets(ns, nelec, &all, &nd);
    int rc = fci_solve_subspace(r, fh, as, all, nd);
    free(all);
    return rc;
}

/* ---------- Sample-Based Quantum Diagonalization (SQD / QSCI) ---------- */

static uint64_t splitmix64(uint64_t *state)
{
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static double rand_u01(uint64_t *state)
{
    return (double)(splitmix64(state) >> 11) * 0x1.0p-53;
}

/* Count alpha (even bits) and beta (odd bits) electrons in a spin-orbital bitstring */
static void count_ab(uint32_t d, int nspin, int *na, int *nb)
{
    int a = 0, b = 0;
    for (int q = 0; q < nspin; q += 2) {
        if ((d >> q) & 1u) a++;
        if (q + 1 < nspin && ((d >> (q + 1)) & 1u)) b++;
    }
    *na = a;
    *nb = b;
}

/* Self-consistent configuration recovery: project a corrupted bitstring into the
   (n_alpha, n_beta) sector weighted by average orbital occupancies n_occ[q]. */
static uint32_t recover_bitstring(uint32_t raw, int nspin, int target_na, int target_nb,
                                  const double *n_occ, uint64_t *rng)
{
    uint32_t out = raw;
    for (int spin = 0; spin < 2; spin++) {
        int target = (spin == 0) ? target_na : target_nb;
        int cur = 0;
        for (int q = spin; q < nspin; q += 2)
            if ((out >> q) & 1u) cur++;
        while (cur > target) {
            /* Too many electrons in this spin sector: flip 1->0 preferring low occupancy orbitals */
            int best_q = -1;
            double best_score = -1.0;
            for (int q = spin; q < nspin; q += 2) {
                if (!((out >> q) & 1u)) continue;
                double w = (1.0 - n_occ[q] + 0.05) * (0.5 + rand_u01(rng));
                if (w > best_score) { best_score = w; best_q = q; }
            }
            if (best_q < 0) break;
            out &= ~(1u << best_q);
            cur--;
        }
        while (cur < target) {
            /* Too few electrons in this spin sector: flip 0->1 preferring high occupancy orbitals */
            int best_q = -1;
            double best_score = -1.0;
            for (int q = spin; q < nspin; q += 2) {
                if ((out >> q) & 1u) continue;
                double w = (n_occ[q] + 0.05) * (0.5 + rand_u01(rng));
                if (w > best_score) { best_score = w; best_q = q; }
            }
            if (best_q < 0) break;
            out |= (1u << best_q);
            cur++;
        }
    }
    return out;
}

int sqd_from_bitstrings(const FermionHam *fh, const ActiveSpace *as,
                        const uint32_t *raw_shots, int nshots, int sc_iters,
                        SQDResult *res)
{
    memset(res, 0, sizeof *res);
    int ns = fh->nspin;
    int target_na = as->nelec_act / 2;
    int target_nb = as->nelec_act - target_na;
    res->n_shots = nshots;

    /* Initial occupancy prior: Hartree-Fock orbitals occupied */
    double n_occ[64] = {0};
    for (int i = 0; i < target_na && 2 * i < ns; i++) n_occ[2 * i] = 0.9;
    for (int i = 0; i < target_nb && 2 * i + 1 < ns; i++) n_occ[2 * i + 1] = 0.9;
    for (int q = 0; q < ns; q++) if (n_occ[q] == 0.0) n_occ[q] = 0.1;

    uint32_t *sub = qc_calloc((size_t)(nshots + 64), sizeof(uint32_t));
    uint64_t rng = 0x123456789abcdef0ULL;

    /* Always include the HF reference determinant */
    uint32_t hf_det = 0;
    for (int i = 0; i < target_na; i++) hf_det |= (1u << (2 * i));
    for (int i = 0; i < target_nb; i++) hf_det |= (1u << (2 * i + 1));

    int valid_cnt = 0;
    for (int s = 0; s < nshots; s++) {
        int na = 0, nb = 0;
        count_ab(raw_shots[s], ns, &na, &nb);
        if (na == target_na && nb == target_nb) valid_cnt++;
    }
    res->n_valid_sector = valid_cnt;

    int n_iters = sc_iters > 0 ? sc_iters : 3;
    FCIResult sub_fci;
    memset(&sub_fci, 0, sizeof sub_fci);

    for (int it = 0; it < n_iters; it++) {
        int nsub = 0;
        sub[nsub++] = hf_det;
        for (int s = 0; s < nshots; s++) {
            int na = 0, nb = 0;
            count_ab(raw_shots[s], ns, &na, &nb);
            if (na == target_na && nb == target_nb) {
                sub[nsub++] = raw_shots[s];
            } else {
                sub[nsub++] = recover_bitstring(raw_shots[s], ns, target_na, target_nb, n_occ, &rng);
            }
        }
        if (it > 0) fci_free(&sub_fci);
        fci_solve_subspace(&sub_fci, fh, as, sub, nsub);

        /* Update average spin-orbital occupancies <psi_SQD | n_q | psi_SQD> */
        memset(n_occ, 0, sizeof n_occ);
        for (int i = 0; i < sub_fci.ndet; i++) {
            double w = sub_fci.ground_vec[i] * sub_fci.ground_vec[i];
            uint32_t d = sub_fci.dets[i];
            for (int q = 0; q < ns; q++)
                if ((d >> q) & 1u) n_occ[q] += w;
        }
    }

    res->n_recovered = nshots - valid_cnt;
    res->subspace_dim = sub_fci.ndet;
    res->energy = sub_fci.ground;
    res->dets = sub_fci.dets;
    res->coeffs = sub_fci.ground_vec;
    free(sub_fci.eigvals);
    free(sub);
    return 0;
}

int sqd_run_simulated(const FermionHam *fh, const ActiveSpace *as, const StateVec *psi,
                      int nshots, double bitflip_prob, int sc_iters, unsigned seed,
                      SQDResult *res)
{
    int ns = fh->nspin;
    size_t dim = psi->dim;
    double *cdf = qc_calloc(dim, sizeof(double));
    double acc = 0.0;
    for (size_t i = 0; i < dim; i++) {
        double p = creal(psi->a[i]) * creal(psi->a[i]) + cimag(psi->a[i]) * cimag(psi->a[i]);
        acc += p;
        cdf[i] = acc;
    }
    if (acc > 0.0) {
        for (size_t i = 0; i < dim; i++) cdf[i] /= acc;
    }

    uint64_t rng = seed ? (uint64_t)seed : 0xCAFEBABE12345678ULL;
    uint32_t *shots = qc_calloc((size_t)nshots, sizeof(uint32_t));
    for (int s = 0; s < nshots; s++) {
        double u = rand_u01(&rng);
        /* binary search in CDF */
        size_t lo = 0, hi = dim - 1;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (cdf[mid] < u) lo = mid + 1;
            else hi = mid;
        }
        uint32_t bitstr = (uint32_t)lo;
        /* inject depolarizing / readout bit-flip noise per qubit */
        if (bitflip_prob > 0.0) {
            for (int q = 0; q < ns; q++) {
                if (rand_u01(&rng) < bitflip_prob) bitstr ^= (1u << q);
            }
        }
        shots[s] = bitstr;
    }
    free(cdf);

    int rc = sqd_from_bitstrings(fh, as, shots, nshots, sc_iters, res);
    /* Also compute what the raw noisy expectation <psi_noisy|H|psi_noisy> would look like */
    if (bitflip_prob > 0.0) {
        /* Depolarizing / bitflip channel attenuates non-identity Pauli weight w by (1 - 2p)^w */
        double e_noisy = 0.0;
        for (size_t k = 0; k < fh->ham.n; k++) {
            uint32_t supp = fh->ham.t[k].x | fh->ham.t[k].z;
            int w = __builtin_popcount(supp);
            double damp = pow(1.0 - 2.0 * bitflip_prob, (double)w);
            double complex ex = sv_expect_pauli(psi, fh->ham.t[k].x, fh->ham.t[k].z);
            e_noisy += creal(fh->ham.t[k].c * ex) * damp;
        }
        res->noisy_estimator_energy = e_noisy;
    } else {
        res->noisy_estimator_energy = sv_expect_sum(psi, &fh->ham);
    }
    free(shots);
    return rc;
}

void sqd_result_free(SQDResult *r)
{
    free(r->dets);
    free(r->coeffs);
    memset(r, 0, sizeof *r);
}

double fci_pauli_restricted_ground(const FermionHam *fh, const FCIResult *fci)
{
    int ns = fh->nspin;
    int nd = fci->ndet;
    static const double complex IPOW4[4] = {1.0, I, -1.0, -I};
    double *M = qc_calloc((size_t)nd * nd, sizeof(double));
    #pragma omp parallel for
    for (int J = 0; J < nd; J++) {
        uint32_t detJ = fci->dets[J];
        for (size_t k = 0; k < fh->ham.n; k++) {
            uint32_t x = fh->ham.t[k].x, z = fh->ham.t[k].z;
            uint32_t detI = detJ ^ x;
            if (__builtin_popcount(detI) != __builtin_popcount(detJ)) {
                if (cabs(fh->ham.t[k].c) > 1e-9) {
                    /* check if it actually connects to outside sector */
                }
                continue;
            }
            int ii = find_det(fci->dets, nd, detI);
            if (ii < 0) continue;
            double sg = (__builtin_popcount(detJ & z) & 1) ? -1.0 : 1.0;
            double complex val = fh->ham.t[k].c * IPOW4[__builtin_popcount(x & z) & 3] * sg;
            M[(size_t)ii * nd + J] += creal(val);
        }
    }
    double *w = qc_calloc((size_t)nd, sizeof(double));
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', nd, M, nd, w);
    double g = w[0];
    free(w);
    free(M);
    (void)ns;
    return g;
}

void fci_free(FCIResult *r)
{
    free(r->eigvals);
    free(r->ground_vec);
    free(r->dets);
    memset(r, 0, sizeof *r);
}
