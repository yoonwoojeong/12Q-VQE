/*
 * Molecule definition and basis-set construction for libcint.
 *
 * The basis file is the Basis Set Exchange NWChem export (spherical). Each
 * "column" of a shell block is one contracted spherical shell; primitives with
 * zero coefficient in a column are dropped. Primitive normalisation is applied
 * for libcint, and each contracted shell is then renormalised using libcint's
 * own overlap so that the diagonal of S is exactly 1.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cint.h>
#include <cint_funcs.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PTR_ENV_OFFSET 20
#define MAX_SHELL_TOTAL 256

void mol_h2_setup(Molecule *m, double bond_angstrom, int charge, int mult)
{
    memset(m, 0, sizeof *m);
    double r = 0.5 * bond_angstrom * QC_BOHR_PER_ANG;
    m->natom = 2;
    m->Z[0] = 1;
    m->Z[1] = 1;
    m->xyz[0][2] = -r;
    m->xyz[1][2] = +r;
    m->charge = charge;
    m->mult = mult;
    m->nelec = 2 - charge;
    double enuc = 0.0;
    for (int i = 0; i < m->natom; i++)
        for (int j = i + 1; j < m->natom; j++) {
            double d = 0;
            for (int k = 0; k < 3; k++) {
                double t = m->xyz[i][k] - m->xyz[j][k];
                d += t * t;
            }
            enuc += (double)(m->Z[i] * m->Z[j]) / sqrt(d);
        }
    m->enuc = enuc;
}

/* ---- parsing of Basis Set Exchange NWChem format ---- */

typedef struct {
    int l;
    int nprim;
    double exp[QC_MAX_SHELL_PRIM];
    double coef[QC_MAX_SHELL_PRIM];
} ParsedShell;

static int l_from_letter(char c)
{
    switch (c) {
    case 'S': return 0;
    case 'P': return 1;
    case 'D': return 2;
    case 'F': return 3;
    case 'G': return 4;
    case 'H': return 5;
    default: return -1;
    }
}

/* Parse all shells for element symbol `sym` from an NWChem BSE file. */
static int parse_nwchem_element(const char *text, const char *sym, ParsedShell *out, int maxout, int *nout)
{
    *nout = 0;
    const char *p = text;
    int in_block = 0;
    int cur_l = -1;
    int ncol = 0;
    /* collect columns: for a block, column j is a list of (exp, coef) */
    static double colexp[16][QC_MAX_SHELL_PRIM];
    static double colcoef[16][QC_MAX_SHELL_PRIM];
    static int colnp[16];
    int nprim_seen = 0;

    while (*p) {
        const char *eol = strchr(p, '\n');
        size_t len = eol ? (size_t)(eol - p) : strlen(p);
        char line[512];
        if (len >= sizeof line) len = sizeof line - 1;
        memcpy(line, p, len);
        line[len] = 0;
        p = eol ? eol + 1 : p + len;

        char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == 0 || *s == '#') continue;

        char tok0[32] = {0}, tok1[32] = {0};
        int nt = sscanf(s, "%31s %31s", tok0, tok1);
        if (strcmp(tok0, "END") == 0) break;
        if (strcmp(tok0, "BASIS") == 0) continue;

        if (nt >= 2 && isalpha((unsigned char)tok0[0]) && strlen(tok1) == 1 && isalpha((unsigned char)tok1[0])) {
            /* new shell block header: "H    S" */
            /* flush previous block */
            if (in_block && strcmp(tok0, sym) == 0) { /* same element, new l */ }
            if (in_block && cur_l >= 0) {
                for (int c = 0; c < ncol; c++) {
                    if (*nout >= maxout) return -1;
                    ParsedShell *ps = &out[(*nout)++];
                    ps->l = cur_l;
                    ps->nprim = 0;
                    for (int k = 0; k < colnp[c]; k++) {
                        if (colcoef[c][k] == 0.0) continue;
                        ps->exp[ps->nprim] = colexp[c][k];
                        ps->coef[ps->nprim] = colcoef[c][k];
                        ps->nprim++;
                    }
                }
            }
            if (in_block && strcmp(tok0, sym) != 0) {
                /* block of a different element: we are done with ours */
            }
            in_block = (strcmp(tok0, sym) == 0);
            cur_l = in_block ? l_from_letter(tok1[0]) : -1;
            ncol = 0;
            nprim_seen = 0;
            memset(colnp, 0, sizeof colnp);
            if (in_block && cur_l < 0) return -2;
            continue;
        }

        if (!in_block) continue;

        /* primitive line: exponent followed by one coefficient per column */
        char *tk = s;
        char *endp = NULL;
        double ex = strtod(tk, &endp);
        if (endp == tk) continue;
        tk = endp;
        int col = 0;
        int nc = 0;
        while (1) {
            while (*tk == ' ' || *tk == '\t') tk++;
            if (*tk == 0 || *tk == '\r' || *tk == '\n') break;
            double cf = strtod(tk, &endp);
            if (endp == tk) break;
            tk = endp;
            if (nprim_seen == 0) ncol++;   /* columns are fixed by the first primitive line */
            if (col >= 16) return -3;
            colexp[col][colnp[col]] = ex;
            colcoef[col][colnp[col]] = cf;
            colnp[col]++;
            col++;
            nc++;
        }
        nprim_seen++;
    }
    /* flush last block */
    if (in_block && cur_l >= 0) {
        for (int c = 0; c < ncol; c++) {
            if (*nout >= maxout) return -1;
            ParsedShell *ps = &out[(*nout)++];
            ps->l = cur_l;
            ps->nprim = 0;
            for (int k = 0; k < colnp[c]; k++) {
                if (colcoef[c][k] == 0.0) continue;
                ps->exp[ps->nprim] = colexp[c][k];
                ps->coef[ps->nprim] = colcoef[c][k];
                ps->nprim++;
            }
        }
    }
    return *nout > 0 ? 0 : -4;
}

/* ---- construction of libcint arrays ---- */

int basis_build(BasisSet *bs, const Molecule *mol, const char *nwchem_path, const char *name)
{
    memset(bs, 0, sizeof *bs);
    snprintf(bs->name, sizeof bs->name, "%s", name);

    size_t len = 0;
    char *text = qc_read_file(nwchem_path, &len);
    if (!text) qc_die("cannot read basis file '%s'", nwchem_path);

    /* Parse shells for every distinct element present (only H is needed here). */
    ParsedShell *elem_shells = qc_calloc(MAX_SHELL_TOTAL, sizeof(ParsedShell));
    int nel_shells = 0;
    int ok = 0;
    /* parse the element used by the molecule; only H is supported by this build */
    for (int a = 0; a < mol->natom; a++) {
        if (mol->Z[a] != 1) qc_die("basis '%s': element Z=%d not available in this build (H only)", name, mol->Z[a]);
    }
    ok = parse_nwchem_element(text, "H", elem_shells, MAX_SHELL_TOTAL, &nel_shells);
    free(text);
    if (ok != 0) qc_die("basis '%s': failed to parse H shells (code %d)", nwchem_path, ok);

    int nshell = nel_shells * mol->natom;
    int nprim_total = 0;
    for (int i = 0; i < nel_shells; i++) {
        for (int a = 0; a < mol->natom; a++) nprim_total += elem_shells[i].nprim;
    }
    (void)nprim_total;

    /* count primitives for env allocation */
    int env_cap = PTR_ENV_OFFSET + 3 * mol->natom + 4 * nshell * QC_MAX_SHELL_PRIM + 16;
    bs->env = qc_calloc((size_t)env_cap, sizeof(double));
    bs->atm = qc_calloc((size_t)6 * mol->natom, sizeof(int));
    bs->bas = qc_calloc((size_t)8 * nshell, sizeof(int));
    bs->ao_loc = qc_calloc((size_t)nshell + 1, sizeof(int));
    bs->shell_atom = qc_calloc((size_t)nshell, sizeof(int));
    bs->shell_l = qc_calloc((size_t)nshell, sizeof(int));
    bs->natm = mol->natom;
    bs->nbas = nshell;

    /* env header */
    bs->env[0] = 50.0;     /* PTR_EXPCUTOFF: prescreening cutoff, ln(threshold) */
    bs->env[1] = 0.0;      /* PTR_COMMON_ORIG */
    bs->env[2] = 0.0;
    bs->env[3] = 0.0;

    int nenv = PTR_ENV_OFFSET;
    /* atoms: coordinates into env */
    for (int a = 0; a < mol->natom; a++) {
        int *atm = bs->atm + 6 * a;
        atm[0] = mol->Z[a];                 /* CHARGE_OF */
        atm[1] = nenv;                      /* PTR_COORD */
        atm[2] = 0;                         /* NUC_MOD_OF: point nucleus */
        atm[3] = 0; atm[4] = 0; atm[5] = 0;
        for (int k = 0; k < 3; k++) bs->env[nenv + k] = mol->xyz[a][k];
        nenv += 3;
    }

    int ish = 0;
    int ao = 0;
    for (int a = 0; a < mol->natom; a++) {
        for (int s = 0; s < nel_shells; s++) {
            const ParsedShell *ps = &elem_shells[s];
            int *bas = bs->bas + 8 * ish;
            int l = ps->l;
            int np = ps->nprim;
            bas[0] = a;                 /* ATOM_OF */
            bas[1] = l;                 /* ANG_OF */
            bas[2] = np;                /* NPRIM_OF */
            bas[3] = 1;                 /* NCTR_OF */
            bas[4] = 0;                 /* KAPPA_OF (spherical is set via cint function choice) */
            bas[5] = nenv;              /* PTR_EXP */
            for (int k = 0; k < np; k++) bs->env[nenv + k] = ps->exp[k];
            nenv += np;
            bas[6] = nenv;              /* PTR_COEFF */
            /* primitive normalisation for libcint: c_k * N_k(l, alpha_k) */
            for (int k = 0; k < np; k++) bs->env[nenv + k] = ps->coef[k] * CINTgto_norm(l, ps->exp[k]);
            nenv += np;
            bas[7] = 0;
            bs->shell_atom[ish] = a;
            bs->shell_l[ish] = l;
            bs->ao_loc[ish] = ao;
            ao += 2 * l + 1;
            ish++;
        }
    }
    bs->ao_loc[nshell] = ao;
    bs->nao = ao;
    bs->nenv = nenv;

    /* renormalise each contracted shell with libcint overlap of its own shell pair */
    for (int i = 0; i < nshell; i++) {
        int l = bs->shell_l[i];
        int nc = 2 * l + 1;
        double *buf = qc_calloc((size_t)nc * nc, sizeof(double));
        int shls[2] = {i, i};
        int dims[2] = {nc, nc};
        int c = 0;
        (void)c;
        int ret = int1e_ovlp_sph(buf, dims, shls, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
        (void)ret;
        double sii = buf[0];
        if (!(sii > 0.0)) qc_die("basis: non-positive self-overlap on shell %d", i);
        double f = 1.0 / sqrt(sii);
        double *coef = bs->env + bs->bas[8 * i + 6];
        int np = bs->bas[8 * i + 2];
        for (int k = 0; k < np; k++) coef[k] *= f;
        free(buf);
    }

    free(elem_shells);
    return 0;
}

void basis_free(BasisSet *bs)
{
    free(bs->atm);
    free(bs->bas);
    free(bs->env);
    free(bs->ao_loc);
    free(bs->shell_atom);
    free(bs->shell_l);
    memset(bs, 0, sizeof *bs);
}
