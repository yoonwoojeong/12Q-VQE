/*
 * Pauli-sum algebra.
 *
 * A Pauli term is  c * P(x,z),  with the canonical Hermitian Pauli
 *      P(x,z) = i^{popcount(x & z)} X^x Z^z
 * so that each qubit carries I (x=0,z=0), X (1,0), Z (0,1) or Y (1,1).
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void pauli_sum_init(PauliSum *ps, int nq)
{
    ps->nq = nq;
    ps->n = 0;
    ps->cap = 64;
    ps->t = qc_calloc(ps->cap, sizeof(PauliTerm));
}

void pauli_sum_free(PauliSum *ps)
{
    free(ps->t);
    ps->t = NULL;
    ps->n = ps->cap = 0;
}

void pauli_sum_add(PauliSum *ps, uint32_t x, uint32_t z, double complex c)
{
    if (ps->n == ps->cap) {
        ps->cap = ps->cap ? 2 * ps->cap : 64;
        ps->t = qc_realloc(ps->t, ps->cap * sizeof(PauliTerm));
    }
    ps->t[ps->n].x = x;
    ps->t[ps->n].z = z;
    ps->t[ps->n].c = c;
    ps->n++;
}

static int cmp_term(const void *a, const void *b)
{
    const PauliTerm *p = a, *q = b;
    if (p->x != q->x) return p->x < q->x ? -1 : 1;
    if (p->z != q->z) return p->z < q->z ? -1 : 1;
    return 0;
}

void pauli_sum_merge(PauliSum *ps, double tol)
{
    if (ps->n == 0) return;
    qsort(ps->t, ps->n, sizeof(PauliTerm), cmp_term);
    size_t w = 0;
    for (size_t r = 0; r < ps->n;) {
        size_t r2 = r;
        double complex c = 0.0;
        while (r2 < ps->n && ps->t[r2].x == ps->t[r].x && ps->t[r2].z == ps->t[r].z) {
            c += ps->t[r2].c;
            r2++;
        }
        if (cabs(c) > tol) {
            ps->t[w].x = ps->t[r].x;
            ps->t[w].z = ps->t[r].z;
            ps->t[w].c = c;
            w++;
        }
        r = r2;
    }
    ps->n = w;
}

void pauli_sum_scale(PauliSum *ps, double complex s)
{
    for (size_t k = 0; k < ps->n; k++) ps->t[k].c *= s;
}

void pauli_sum_append(PauliSum *dst, const PauliSum *src, double complex s)
{
    for (size_t k = 0; k < src->n; k++) pauli_sum_add(dst, src->t[k].x, src->t[k].z, s * src->t[k].c);
}

/* phase factors i^k, k mod 4 */
static const double complex IPOW[4] = {1.0, I, -1.0, -I};

PauliSum pauli_mul(const PauliSum *a, const PauliSum *b, double tol)
{
    PauliSum out;
    pauli_sum_init(&out, a->nq);
    for (size_t i = 0; i < a->n; i++) {
        const PauliTerm *ta = &a->t[i];
        int ka = __builtin_popcount(ta->x & ta->z);
        for (size_t j = 0; j < b->n; j++) {
            const PauliTerm *tb = &b->t[j];
            int kb = __builtin_popcount(tb->x & tb->z);
            uint32_t x = ta->x ^ tb->x;
            uint32_t z = ta->z ^ tb->z;
            int kc = __builtin_popcount(x & z);
            /* X^xa Z^za X^xb Z^zb = (-1)^{|za & xb|} X^x Z^z */
            int sign = (__builtin_popcount(ta->z & tb->x) & 1) ? -1 : 1;
            int e = ((ka + kb - kc) % 4 + 4) % 4;
            double complex c = ta->c * tb->c * IPOW[e] * (double)sign;
            if (c != 0.0) pauli_sum_add(&out, x, z, c);
        }
    }
    pauli_sum_merge(&out, tol);
    return out;
}

int pauli_label(uint32_t x, uint32_t z, int nq, char *out)
{
    for (int q = 0; q < nq; q++) {
        int xb = (x >> q) & 1, zb = (z >> q) & 1;
        out[q] = (xb && zb) ? 'Y' : (xb ? 'X' : (zb ? 'Z' : 'I'));
    }
    out[nq] = 0;
    return nq;
}

int pauli_is_hermitian(const PauliSum *ps, double tol, double *max_imag)
{
    double mi = 0.0;
    for (size_t k = 0; k < ps->n; k++) {
        double im = fabs(cimag(ps->t[k].c));
        if (im > mi) mi = im;
    }
    if (max_imag) *max_imag = mi;
    return mi <= tol;
}
