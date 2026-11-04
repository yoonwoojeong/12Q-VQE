/*
 * Molecule definition and basis-set construction for libcint.
 *
 * The basis file is the Basis Set Exchange NWChem export (spherical). Each
 * "column" of a shell block is one contracted spherical shell; primitives with
 * zero coefficient in a column are dropped. SP blocks split into an S and a P
 * shell. Primitive normalisation is applied for libcint, and each contracted
 * shell is then renormalised using libcint's own overlap so that the diagonal
 * of S is exactly 1.
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
#define MAX_SHELL_TOTAL 512

static void compute_enuc(Molecule *m)
{
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
    compute_enuc(m);
}

int mol_setup_by_name(Molecule *m, const char *name, double bond_angstrom)
{
    memset(m, 0, sizeof *m);
    double r = bond_angstrom * QC_BOHR_PER_ANG;
    if (strcmp(name, "H2") == 0 || strcmp(name, "h2") == 0) {
        mol_h2_setup(m, bond_angstrom, 0, 1);
        return 0;
    }
    if (strcmp(name, "HeH+") == 0 || strcmp(name, "heh+") == 0) {
        m->natom = 2;
        m->Z[0] = 2; m->Z[1] = 1;
        m->xyz[0][2] = -0.5 * r;
        m->xyz[1][2] = +0.5 * r;
        m->charge = 1; m->mult = 1; m->nelec = 2;
        compute_enuc(m);
        return 0;
    }
    if (strcmp(name, "LiH") == 0 || strcmp(name, "lih") == 0) {
        m->natom = 2;
        m->Z[0] = 3; m->Z[1] = 1;
        m->xyz[0][2] = -0.5 * r;
        m->xyz[1][2] = +0.5 * r;
        m->charge = 0; m->mult = 1; m->nelec = 4;
        compute_enuc(m);
        return 0;
    }
    if (strcmp(name, "BeH2") == 0 || strcmp(name, "beh2") == 0) {
        m->natom = 3;
        m->Z[0] = 4; m->Z[1] = 1; m->Z[2] = 1;
        m->xyz[0][2] = 0.0;
        m->xyz[1][2] = -r;
        m->xyz[2][2] = +r;
        m->charge = 0; m->mult = 1; m->nelec = 6;
        compute_enuc(m);
        return 0;
    }
    if (strcmp(name, "H4") == 0 || strcmp(name, "h4") == 0) {
        /* Linear H4 chain with nearest-neighbour spacing bond_angstrom */
        m->natom = 4;
        for (int i = 0; i < 4; i++) {
            m->Z[i] = 1;
            m->xyz[i][2] = ((double)i - 1.5) * r;
        }
        m->charge = 0; m->mult = 1; m->nelec = 4;
        compute_enuc(m);
        return 0;
    }
    return -1;
}

static const char *ELEM_SYMBOLS[] = {
    NULL,
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr"
};

const char *elem_symbol(int Z)
{
    if (Z >= 1 && Z <= 36) return ELEM_SYMBOLS[Z];
    return NULL;
}

int elem_z_from_symbol(const char *sym)
{
    if (!sym || !*sym) return -1;
    char norm[4] = {0};
    norm[0] = (char)toupper((unsigned char)sym[0]);
    if (sym[1] && isalpha((unsigned char)sym[1])) norm[1] = (char)tolower((unsigned char)sym[1]);
    for (int z = 1; z <= 36; z++) {
        if (strcmp(norm, ELEM_SYMBOLS[z]) == 0) return z;
    }
    int znum = atoi(sym);
    if (znum >= 1 && znum <= 36) return znum;
    return -1;
}

int mol_load_xyz(Molecule *m, const char *xyz_path, int charge, int mult)
{
    size_t len = 0;
    char *txt = qc_read_file(xyz_path, &len);
    if (!txt) return -1;
    memset(m, 0, sizeof *m);
    const char *p = txt;
    /* line 1: natom */
    const char *eol = strchr(p, '\n');
    int nat = atoi(p);
    if (nat <= 0 || nat > QC_MAX_ATOM || !eol) { free(txt); return -2; }
    p = eol + 1;
    /* line 2: comment (may optionally contain charge=<int> mult=<int>) */
    eol = strchr(p, '\n');
    if (!eol) { free(txt); return -3; }
    p = eol + 1;

    m->natom = nat;
    int sum_z = 0;
    for (int a = 0; a < nat; a++) {
        while (*p == '\r' || *p == '\n') p++;
        char sym[16] = {0};
        double x = 0, y = 0, z = 0;
        int nread = sscanf(p, "%15s %lf %lf %lf", sym, &x, &y, &z);
        if (nread != 4) { free(txt); return -4; }
        int Z = elem_z_from_symbol(sym);
        if (Z <= 0) { free(txt); return -5; }
        m->Z[a] = Z;
        m->xyz[a][0] = x * QC_BOHR_PER_ANG;
        m->xyz[a][1] = y * QC_BOHR_PER_ANG;
        m->xyz[a][2] = z * QC_BOHR_PER_ANG;
        sum_z += Z;
        eol = strchr(p, '\n');
        p = eol ? eol + 1 : p + strlen(p);
    }
    free(txt);
    m->charge = charge;
    m->mult = mult > 0 ? mult : 1;
    m->nelec = sum_z - charge;
    compute_enuc(m);
    return 0;
}

/* ---- parsing of Basis Set Exchange NWChem format ---- */

typedef struct {
    int l;
    int nprim;
    double exp[QC_MAX_SHELL_PRIM];
    double coef[QC_MAX_SHELL_PRIM];
} ParsedShell;

static int l_from_token(const char *tok)
{
    if (strcmp(tok, "SP") == 0) return -10; /* special SP marker: col 0 = S, col 1 = P */
    if (strlen(tok) != 1) return -1;
    switch (tok[0]) {
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

        if (nt >= 2 && isalpha((unsigned char)tok0[0]) && isalpha((unsigned char)tok1[0])) {
            /* flush previous block if it belonged to our element */
            if (in_block && (cur_l >= 0 || cur_l == -10)) {
                for (int c = 0; c < ncol; c++) {
                    if (*nout >= maxout) return -1;
                    ParsedShell *ps = &out[(*nout)++];
                    ps->l = (cur_l == -10) ? c : cur_l;
                    ps->nprim = 0;
                    for (int k = 0; k < colnp[c]; k++) {
                        if (colcoef[c][k] == 0.0) continue;
                        ps->exp[ps->nprim] = colexp[c][k];
                        ps->coef[ps->nprim] = colcoef[c][k];
                        ps->nprim++;
                    }
                }
            }
            in_block = (strcmp(tok0, sym) == 0);
            cur_l = in_block ? l_from_token(tok1) : -1;
            ncol = 0;
            nprim_seen = 0;
            memset(colnp, 0, sizeof colnp);
            if (in_block && cur_l < 0 && cur_l != -10) return -2;
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
        }
        nprim_seen++;
    }
    /* flush last block */
    if (in_block && (cur_l >= 0 || cur_l == -10)) {
        for (int c = 0; c < ncol; c++) {
            if (*nout >= maxout) return -1;
            ParsedShell *ps = &out[(*nout)++];
            ps->l = (cur_l == -10) ? c : cur_l;
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

    /* Parse shells for each atom */
    ParsedShell *atom_shells[QC_MAX_ATOM];
    int natom_shells[QC_MAX_ATOM];
    int nshell = 0;
    for (int a = 0; a < mol->natom; a++) {
        const char *sym = elem_symbol(mol->Z[a]);
        if (!sym) qc_die("basis '%s': unsupported atomic number Z=%d", name, mol->Z[a]);
        atom_shells[a] = qc_calloc(MAX_SHELL_TOTAL, sizeof(ParsedShell));
        int ok = parse_nwchem_element(text, sym, atom_shells[a], MAX_SHELL_TOTAL, &natom_shells[a]);
        if (ok != 0) qc_die("basis '%s': failed to parse %s shells (code %d)", nwchem_path, sym, ok);
        nshell += natom_shells[a];
    }
    free(text);

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
        for (int s = 0; s < natom_shells[a]; s++) {
            const ParsedShell *ps = &atom_shells[a][s];
            int *bas = bs->bas + 8 * ish;
            int l = ps->l;
            int np = ps->nprim;
            bas[0] = a;                 /* ATOM_OF */
            bas[1] = l;                 /* ANG_OF */
            bas[2] = np;                /* NPRIM_OF */
            bas[3] = 1;                 /* NCTR_OF */
            bas[4] = 0;                 /* KAPPA_OF */
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
        free(atom_shells[a]);
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
        (void)int1e_ovlp_sph(buf, dims, shls, bs->atm, bs->natm, bs->bas, bs->nbas, bs->env, NULL, NULL);
        double sii = buf[0];
        if (!(sii > 0.0)) qc_die("basis: non-positive self-overlap on shell %d", i);
        double f = 1.0 / sqrt(sii);
        double *coef = bs->env + bs->bas[8 * i + 6];
        int np = bs->bas[8 * i + 2];
        for (int k = 0; k < np; k++) coef[k] *= f;
        free(buf);
    }

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
