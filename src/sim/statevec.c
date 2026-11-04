/*
 * Dense statevector simulator (up to 20 qubits, OpenMP-accelerated).
 * Basis index bit k = occupation of qubit k.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static const double complex IPOW4[4] = {1.0, I, -1.0, -I};

int sv_init(StateVec *sv, int nq)
{
    if (nq < 1 || nq > QC_MAX_QUBIT_SV) qc_die("statevector: %d qubits outside 1..%d", nq, QC_MAX_QUBIT_SV);
    sv->nq = nq;
    sv->dim = (size_t)1 << nq;
    sv->a = qc_calloc(sv->dim, sizeof(double complex));
    sv->a[0] = 1.0;
    return 0;
}

void sv_free(StateVec *sv)
{
    free(sv->a);
    sv->a = NULL;
    sv->dim = 0;
}

void sv_set_basis(StateVec *sv, size_t index)
{
    memset(sv->a, 0, sizeof(double complex) * sv->dim);
    sv->a[index] = 1.0;
}

void sv_copy(StateVec *dst, const StateVec *src)
{
    if (dst->dim != src->dim) qc_die("sv_copy: dimension mismatch");
    memcpy(dst->a, src->a, sizeof(double complex) * src->dim);
}

/* psi <- exp(i phi P) psi, P = i^k X^x Z^z Hermitian, P^2 = 1 */
void sv_apply_pauli_rotation(StateVec *sv, uint32_t x, uint32_t z, double phi)
{
    double complex c = cos(phi);
    double complex s = I * sin(phi);
    int k = __builtin_popcount(x & z) & 3;
    double complex ik = IPOW4[k];
    size_t dim = sv->dim;
    double complex *a = sv->a;
    if (x == 0) {
        /* diagonal: P psi_i = i^k (-1)^{pc(i&z)} psi_i */
        #pragma omp parallel for if(dim >= 16384)
        for (size_t i = 0; i < dim; i++) {
            double sg = (__builtin_popcount((unsigned)(i & z)) & 1) ? -1.0 : 1.0;
            double complex p = ik * sg;
            a[i] = c * a[i] + s * p * a[i];
        }
        return;
    }
    #pragma omp parallel for if(dim >= 16384)
    for (size_t i = 0; i < dim; i++) {
        size_t j = i ^ x;
        if (j < i) continue;          /* handle each pair once */
        double complex ai = a[i], aj = a[j];
        double sgj = (__builtin_popcount((unsigned)(j & z)) & 1) ? -1.0 : 1.0;   /* (P psi)_i uses psi_j */
        double sgi = (__builtin_popcount((unsigned)(i & z)) & 1) ? -1.0 : 1.0;   /* (P psi)_j uses psi_i */
        a[i] = c * ai + s * ik * sgj * aj;
        a[j] = c * aj + s * ik * sgi * ai;
    }
}

/* out <- (sum_t c_t P_t) in */
void sv_apply_sum(const StateVec *in, const PauliSum *ps, StateVec *out)
{
    size_t dim = in->dim;
    #pragma omp parallel for if(dim >= 4096)
    for (size_t i = 0; i < dim; i++) {
        double complex acc = 0.0;
        for (size_t k = 0; k < ps->n; k++) {
            uint32_t x = ps->t[k].x, z = ps->t[k].z;
            double complex ph = ps->t[k].c * IPOW4[__builtin_popcount(x & z) & 3];
            size_t j = i ^ x;
            double sg = (__builtin_popcount((unsigned)(j & z)) & 1) ? -1.0 : 1.0;
            acc += ph * sg * in->a[j];
        }
        out->a[i] = acc;
    }
}

void sv_apply_1q(StateVec *sv, int q, const double complex u[2][2])
{
    size_t bit = (size_t)1 << q;
    size_t dim = sv->dim;
    double complex u00 = u[0][0], u01 = u[0][1], u10 = u[1][0], u11 = u[1][1];
    #pragma omp parallel for if(dim >= 16384)
    for (size_t i = 0; i < dim; i++) {
        if (i & bit) continue;
        size_t j = i | bit;
        double complex a = sv->a[i], b = sv->a[j];
        sv->a[i] = u00 * a + u01 * b;
        sv->a[j] = u10 * a + u11 * b;
    }
}

void sv_apply_cx(StateVec *sv, int c, int t)
{
    size_t cb = (size_t)1 << c, tb = (size_t)1 << t;
    size_t dim = sv->dim;
    #pragma omp parallel for if(dim >= 16384)
    for (size_t i = 0; i < dim; i++) {
        if ((i & cb) && !(i & tb)) {
            double complex tmp = sv->a[i];
            sv->a[i] = sv->a[i | tb];
            sv->a[i | tb] = tmp;
        }
    }
}

void sv_apply_cz(StateVec *sv, int c, int t)
{
    size_t cb = (size_t)1 << c, tb = (size_t)1 << t;
    size_t dim = sv->dim;
    #pragma omp parallel for if(dim >= 16384)
    for (size_t i = 0; i < dim; i++)
        if ((i & cb) && (i & tb)) sv->a[i] = -sv->a[i];
}

double complex sv_expect_pauli(const StateVec *sv, uint32_t x, uint32_t z)
{
    int k = __builtin_popcount(x & z) & 3;
    double complex ik = IPOW4[k];
    double sr = 0.0, si = 0.0;
    size_t dim = sv->dim;
    #pragma omp parallel for reduction(+:sr,si) if(dim >= 16384)
    for (size_t i = 0; i < dim; i++) {
        size_t j = i ^ x;
        double sg = (__builtin_popcount((unsigned)(j & z)) & 1) ? -1.0 : 1.0;
        double complex val = conj(sv->a[i]) * ik * sg * sv->a[j];
        sr += creal(val);
        si += cimag(val);
    }
    return sr + I * si;
}

double sv_expect_sum(const StateVec *sv, const PauliSum *ps)
{
    double e = 0.0;
    size_t nt = ps->n;
    #pragma omp parallel for reduction(+:e) if(nt >= 16)
    for (size_t t = 0; t < nt; t++) {
        uint32_t x = ps->t[t].x, z = ps->t[t].z;
        int k = __builtin_popcount(x & z) & 3;
        double complex ik = IPOW4[k];
        double complex s = 0.0;
        for (size_t i = 0; i < sv->dim; i++) {
            size_t j = i ^ x;
            double sg = (__builtin_popcount((unsigned)(j & z)) & 1) ? -1.0 : 1.0;
            s += conj(sv->a[i]) * ik * sg * sv->a[j];
        }
        e += creal(ps->t[t].c * s);
    }
    return e;
}

double sv_norm(const StateVec *sv)
{
    double s = 0.0;
    for (size_t i = 0; i < sv->dim; i++) s += creal(sv->a[i] * conj(sv->a[i]));
    return sqrt(s);
}
