/*
 * Active-space construction: AO -> MO integral transformation for a window of
 * spatial orbitals [first, first+norb). Doubly occupied orbitals below the
 * window are inactive and are folded into an effective one-body operator and a
 * constant energy (E_core).
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cblas.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* out[a][b][c][d] = sum (mu nu|la si) Ca[mu][a] Cb[nu][b] Cc[la][c] Cd[si][d] (chemist) */
static double *ao2mo4(int n, const double *eri, const double *Ca, int na, const double *Cb, int nb,
                      const double *Cc, int nc, const double *Cd, int nd)
{
    size_t n2 = (size_t)n * n;
    size_t n3 = n2 * (size_t)n;
    /* step 1: contract first index mu -> a: T1[a, nu*la*si] = Ca^T[a, mu] * eri[mu, nu*la*si] */
    double *T1 = qc_calloc((size_t)na * n3, sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                na, (int)n3, n, 1.0, Ca, na, eri, (int)n3, 0.0, T1, (int)n3);

    /* step 2: second index nu -> b: T2[a, b, la*si] = sum_nu Cb^T[b, nu] * T1[a, nu, la*si] */
    double *T2 = qc_calloc((size_t)na * nb * n2, sizeof(double));
    for (int a = 0; a < na; a++) {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                    nb, (int)n2, n, 1.0, Cb, nb, T1 + (size_t)a * n3, (int)n2,
                    0.0, T2 + (size_t)a * nb * n2, (int)n2);
    }
    free(T1);

    /* step 3: la -> c : T3[a*nb + b, c, si] = sum_la Cc^T[c, la] * T2[a*nb + b, la, si] */
    double *T3 = qc_calloc((size_t)na * nb * nc * n, sizeof(double));
    for (int ab = 0; ab < na * nb; ab++) {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                    nc, n, n, 1.0, Cc, nc, T2 + (size_t)ab * n2, n,
                    0.0, T3 + (size_t)ab * nc * n, n);
    }
    free(T2);

    /* step 4: si -> d : out[a*nb*nc + c, d] = sum_si T3[a*nb*nc + c, si] * Cd[si, d] */
    double *out = qc_calloc((size_t)na * nb * nc * nd, sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                na * nb * nc, nd, n, 1.0, T3, n, Cd, nd, 0.0, out, nd);
    free(T3);
    return out;
}

/* one-electron MO matrix: out[p][q] = sum Ca[mu][p] Cb[nu][q] M[mu][nu] */
static void mo1e(int n, const double *M, const double *Ca, int na, const double *Cb, int nb, double *out)
{
    double *tmp = qc_calloc((size_t)n * nb, sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, nb, n, 1.0, M, n, Cb, nb, 0.0, tmp, nb);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, na, nb, n, 1.0, Ca, na, tmp, nb, 0.0, out, nb);
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

int cas_system_load(CasSystem *cs, int nel, int norb, double bond)
{
    memset(cs, 0, sizeof *cs);
    char ddir[512], pth[1024];
    qc_data_dir(ddir, sizeof ddir);
    snprintf(pth, sizeof pth, "%s/basis/cc-pVQZ.nw", ddir);
    mol_h2_setup(&cs->mol, bond, 0, 1);
    BasisSet bs;
    basis_build(&bs, &cs->mol, pth, "cc-pVQZ");
    AOInts ints;
    char cdir[512], cache[1024];
    qc_cache_dir(cdir, sizeof cdir);
    qc_mkdir_p(cdir);
    snprintf(cache, sizeof cache, "%s/ao_cc-pVQZ_H2_%.4f.bin", cdir, bond);
    if (ao_integrals_read(&ints, cache, bs.nao) != 0) {
        ao_integrals_compute(&ints, &bs, &cs->mol);
        ao_integrals_write(&ints, cache, "cc-pVQZ");
    }
    rhf_run(&cs->rhf, &ints, &cs->mol, 200, 1e-10, 1e-8, 0);
    int first = (cs->mol.nelec - nel) / 2;
    active_space_build(&cs->as, &ints, &cs->rhf, &cs->mol, nel, norb, first);
    fermion_build_hamiltonian(&cs->fh, &cs->as);
    ao_integrals_free(&ints);
    basis_free(&bs);
    return 0;
}

void cas_system_free(CasSystem *cs)
{
    fermion_free(&cs->fh);
    active_space_free(&cs->as);
    rhf_free(&cs->rhf);
    memset(cs, 0, sizeof *cs);
}

int fcidump_write(const ActiveSpace *as, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    int n = as->nact;
    fprintf(f, " &FCI NORB=%d,NELEC=%d,MS2=0,\n", n, as->nelec_act);
    fprintf(f, "  ORBSYM=");
    for (int i = 0; i < n; i++) fprintf(f, "1,");
    fprintf(f, "\n  ISYM=1,\n &END\n");

    double tol = 1e-15;
    /* 2-electron integrals (p q | r s) with 1-based canonical indices p>=q, r>=s, pq>=rs */
    for (int p = 0; p < n; p++)
        for (int q = 0; q <= p; q++) {
            int pq = p * (p + 1) / 2 + q;
            for (int r = 0; r < n; r++)
                for (int s = 0; s <= r; s++) {
                    int rs = r * (r + 1) / 2 + s;
                    if (pq < rs) continue;
                    double v = as->g[(((size_t)p * n + q) * n + r) * n + s];
                    if (fabs(v) > tol)
                        fprintf(f, " %23.16e %3d %3d %3d %3d\n", v, p + 1, q + 1, r + 1, s + 1);
                }
        }
    /* 1-electron integrals h_pq with p>=q */
    for (int p = 0; p < n; p++)
        for (int q = 0; q <= p; q++) {
            double v = as->h[(size_t)p * n + q];
            if (fabs(v) > tol)
                fprintf(f, " %23.16e %3d %3d   0   0\n", v, p + 1, q + 1);
        }
    /* core energy */
    fprintf(f, " %23.16e   0   0   0   0\n", as->ecore);
    fclose(f);
    return 0;
}

static int parse_int_after(const char *s, const char *key, int dflt)
{
    const char *p = strstr(s, key);
    if (!p) return dflt;
    p += strlen(key);
    while (*p == ' ' || *p == '=') p++;
    return atoi(p);
}

int fcidump_read(ActiveSpace *as, const char *path)
{
    size_t len = 0;
    char *txt = qc_read_file(path, &len);
    if (!txt) return -1;
    memset(as, 0, sizeof *as);

    /* Find header end (&END or / ) */
    const char *body = strstr(txt, "&END");
    if (body) body += 4;
    else {
        body = strstr(txt, "&end");
        if (body) body += 4;
        else {
            const char *sl = strchr(txt, '/');
            if (sl) body = sl + 1;
        }
    }
    if (!body) { free(txt); return -2; }

    int norb = parse_int_after(txt, "NORB", 0);
    if (!norb) norb = parse_int_after(txt, "norb", 0);
    int nelec = parse_int_after(txt, "NELEC", 0);
    if (!nelec) nelec = parse_int_after(txt, "nelec", 0);
    if (norb <= 0 || nelec < 0 || nelec > 2 * norb) { free(txt); return -3; }

    as->nact = norb;
    as->nao = norb;
    as->nelec_act = nelec;
    as->h = qc_calloc((size_t)norb * norb, sizeof(double));
    as->g = qc_calloc((size_t)norb * norb * norb * norb, sizeof(double));

    const char *p = body;
    while (*p) {
        while (*p == '\r' || *p == '\n' || *p == ' ' || *p == '\t') p++;
        if (!*p) break;
        const char *eol = strchr(p, '\n');
        double val = 0.0;
        int i = 0, j = 0, k = 0, l = 0;
        if (sscanf(p, "%lf %d %d %d %d", &val, &i, &j, &k, &l) == 5) {
            if (i == 0 && j == 0 && k == 0 && l == 0) {
                as->ecore = val;
            } else if (k == 0 && l == 0 && i >= 1 && i <= norb && j >= 1 && j <= norb) {
                int ii = i - 1, jj = j - 1;
                as->h[(size_t)ii * norb + jj] = val;
                as->h[(size_t)jj * norb + ii] = val;
            } else if (i >= 1 && i <= norb && j >= 1 && j <= norb && k >= 1 && k <= norb && l >= 1 && l <= norb) {
                int a = i - 1, b = j - 1, c = k - 1, d = l - 1;
                /* 8-fold chemist symmetry (ab|cd) */
                int perms[8][4] = {
                    {a, b, c, d}, {b, a, c, d}, {a, b, d, c}, {b, a, d, c},
                    {c, d, a, b}, {d, c, a, b}, {c, d, b, a}, {d, c, b, a}
                };
                for (int m = 0; m < 8; m++) {
                    size_t idx = (((size_t)perms[m][0] * norb + perms[m][1]) * norb + perms[m][2]) * norb + perms[m][3];
                    as->g[idx] = val;
                }
            }
        }
        p = eol ? eol + 1 : p + strlen(p);
    }
    free(txt);
    return 0;
}
