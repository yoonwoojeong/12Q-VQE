/*
 * AO integrals through libcint (spherical Gaussians).
 *
 * libcint returns each shell-block in Fortran order (first index fastest).
 * Block (i,j) of a one-electron integral therefore sits at buf[a + ni*b] for
 * AO a in shell i and AO b in shell j.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cint.h>
#include <cint_funcs.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int nao;
    int magic;
} CacheHeader;

#define AOINT_MAGIC 0x41304e54  /* "AO NT" */

static void scatter_1e(double *M, int nao, const BasisSet *bs, int i, int j, const double *buf)
{
    int ni = 2 * bs->shell_l[i] + 1;
    int nj = 2 * bs->shell_l[j] + 1;
    int oi = bs->ao_loc[i];
    int oj = bs->ao_loc[j];
    for (int a = 0; a < ni; a++)
        for (int b = 0; b < nj; b++) {
            double v = buf[a + ni * b];
            M[(size_t)(oi + a) * nao + (oj + b)] = v;
            M[(size_t)(oj + b) * nao + (oi + a)] = v;
        }
}

int ao_integrals_compute(AOInts *ints, const BasisSet *bs, const Molecule *mol)
{
    (void)mol;
    int nao = bs->nao;
    int nsh = bs->nbas;
    memset(ints, 0, sizeof *ints);
    ints->nao = nao;
    ints->S = qc_calloc((size_t)nao * nao, sizeof(double));
    ints->T = qc_calloc((size_t)nao * nao, sizeof(double));
    ints->V = qc_calloc((size_t)nao * nao, sizeof(double));
    ints->H = qc_calloc((size_t)nao * nao, sizeof(double));
    ints->eri = qc_calloc((size_t)nao * nao * nao * nao, sizeof(double));

    int maxbuf = 0;
    for (int i = 0; i < nsh; i++)
        for (int j = 0; j < nsh; j++) {
            int n = (2 * bs->shell_l[i] + 1) * (2 * bs->shell_l[j] + 1);
            if (n > maxbuf) maxbuf = n;
        }
    double *buf = qc_calloc((size_t)maxbuf, sizeof(double));

    for (int i = 0; i < nsh; i++)
        for (int j = 0; j <= i; j++) {
            int shls[2] = {i, j};
            int ni = 2 * bs->shell_l[i] + 1;
            int nj = 2 * bs->shell_l[j] + 1;
            (void)ni; (void)nj;
            int dims[2] = {0, 0};
            (void)dims;
            memset(buf, 0, sizeof(double) * (size_t)maxbuf);
            (void)int1e_ovlp_sph(buf, NULL, shls, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            scatter_1e(ints->S, nao, bs, i, j, buf);
            memset(buf, 0, sizeof(double) * (size_t)maxbuf);
            (void)int1e_kin_sph(buf, NULL, shls, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            scatter_1e(ints->T, nao, bs, i, j, buf);
            memset(buf, 0, sizeof(double) * (size_t)maxbuf);
            (void)int1e_nuc_sph(buf, NULL, shls, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
            scatter_1e(ints->V, nao, bs, i, j, buf);
        }
    free(buf);

    for (size_t k = 0; k < (size_t)nao * nao; k++) ints->H[k] = ints->T[k] + ints->V[k];

    /* two-electron integrals (ij|kl) over shell quartets with 8-fold symmetry */
    CINTOpt *opt = NULL;
    int2e_optimizer(&opt, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env);
    int max4 = 0;
    for (int i = 0; i < nsh; i++)
        for (int j = 0; j < nsh; j++)
            for (int k = 0; k < nsh; k++)
                for (int l = 0; l < nsh; l++) {
                    int n = (2 * bs->shell_l[i] + 1) * (2 * bs->shell_l[j] + 1) * (2 * bs->shell_l[k] + 1) * (2 * bs->shell_l[l] + 1);
                    if (n > max4) max4 = n;
                }
    double *blk = qc_calloc((size_t)max4, sizeof(double));
    size_t N = (size_t)nao;
    double *eri = ints->eri;
    for (int i = 0; i < nsh; i++)
        for (int j = 0; j <= i; j++)
            for (int k = 0; k < nsh; k++)
                for (int l = 0; l <= k; l++) {
                    if (i * (i + 1) / 2 + j < k * (k + 1) / 2 + l) continue;
                    int shls[4] = {i, j, k, l};
                    memset(blk, 0, sizeof(double) * (size_t)max4);
                    (void)int2e_sph(blk, NULL, shls, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, opt, NULL);
                    int ni = 2 * bs->shell_l[i] + 1, nj = 2 * bs->shell_l[j] + 1;
                    int nk = 2 * bs->shell_l[k] + 1, nl = 2 * bs->shell_l[l] + 1;
                    int oi = bs->ao_loc[i], oj = bs->ao_loc[j], ok = bs->ao_loc[k], ol = bs->ao_loc[l];
                    for (int a = 0; a < ni; a++)
                        for (int b = 0; b < nj; b++)
                            for (int c = 0; c < nk; c++)
                                for (int d = 0; d < nl; d++) {
                                    double v = blk[a + ni * (b + nj * (c + nk * d))];
                                    size_t p = oi + a, q = oj + b, rr = ok + c, s = ol + d;
                                    eri[((p * N + q) * N + rr) * N + s] = v;
                                    eri[((q * N + p) * N + rr) * N + s] = v;
                                    eri[((p * N + q) * N + s) * N + rr] = v;
                                    eri[((q * N + p) * N + s) * N + rr] = v;
                                    eri[((rr * N + s) * N + p) * N + q] = v;
                                    eri[((s * N + rr) * N + p) * N + q] = v;
                                    eri[((rr * N + s) * N + q) * N + p] = v;
                                    eri[((s * N + rr) * N + q) * N + p] = v;
                                }
                }
    free(blk);
    CINTdel_2e_optimizer(&opt);
    return 0;
}

void ao_integrals_free(AOInts *ints)
{
    free(ints->S);
    free(ints->T);
    free(ints->V);
    free(ints->H);
    free(ints->eri);
    memset(ints, 0, sizeof *ints);
}

int ao_integrals_write(const AOInts *ints, const char *path, const char *tag)
{
    (void)tag;
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    CacheHeader h = {ints->nao, AOINT_MAGIC};
    size_t n2 = (size_t)ints->nao * ints->nao;
    size_t n4 = n2 * n2;
    int ok = fwrite(&h, sizeof h, 1, f) == 1;
    ok = ok && fwrite(ints->S, sizeof(double), n2, f) == n2;
    ok = ok && fwrite(ints->T, sizeof(double), n2, f) == n2;
    ok = ok && fwrite(ints->V, sizeof(double), n2, f) == n2;
    ok = ok && fwrite(ints->H, sizeof(double), n2, f) == n2;
    ok = ok && fwrite(ints->eri, sizeof(double), n4, f) == n4;
    fclose(f);
    return ok ? 0 : -1;
}

int ao_integrals_read(AOInts *ints, const char *path, int nao_expected)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    CacheHeader h;
    if (fread(&h, sizeof h, 1, f) != 1 || h.magic != AOINT_MAGIC || h.nao != nao_expected) {
        fclose(f);
        return -1;
    }
    int nao = h.nao;
    memset(ints, 0, sizeof *ints);
    ints->nao = nao;
    size_t n2 = (size_t)nao * nao, n4 = n2 * n2;
    ints->S = qc_calloc(n2, sizeof(double));
    ints->T = qc_calloc(n2, sizeof(double));
    ints->V = qc_calloc(n2, sizeof(double));
    ints->H = qc_calloc(n2, sizeof(double));
    ints->eri = qc_calloc(n4, sizeof(double));
    int ok = fread(ints->S, sizeof(double), n2, f) == n2;
    ok = ok && fread(ints->T, sizeof(double), n2, f) == n2;
    ok = ok && fread(ints->V, sizeof(double), n2, f) == n2;
    ok = ok && fread(ints->H, sizeof(double), n2, f) == n2;
    ok = ok && fread(ints->eri, sizeof(double), n4, f) == n4;
    fclose(f);
    if (!ok) { ao_integrals_free(ints); return -1; }
    return 0;
}
