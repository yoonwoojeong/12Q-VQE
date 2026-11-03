/*
 * qchem-vqe command-line interface.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cjson/cJSON.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef QC_DATA_DIR
#define QC_DATA_DIR "data"
#endif
#define CACHE_DIR "build/cache"
#define RESULTS_DIR "results"

typedef struct {
    char basis[64];
    char results_dir[512];
    int electrons;      /* active electrons */
    int orbitals;       /* active spatial orbitals */
    int first;          /* first active MO, -1 = default */
    int maxeval;
    int verbose;
    double bond;        /* Angstrom */
} Opts;

typedef struct {
    Molecule mol;
    BasisSet bs;
    AOInts ints;
    RHFResult rhf;
    ActiveSpace as;
    FermionHam fh;
    int have_ints, have_rhf;
} Ctx;

static void opts_default(Opts *o)
{
    memset(o, 0, sizeof *o);
    snprintf(o->basis, sizeof o->basis, "cc-pVQZ");
    snprintf(o->results_dir, sizeof o->results_dir, "%s", RESULTS_DIR);
    o->electrons = 2;
    o->orbitals = 6;
    o->first = -1;
    o->maxeval = 4000;
    o->bond = 0.735;
}

static const char *opt_value(int argc, char **argv, int *i, const char *name)
{
    if (strcmp(argv[*i], name) != 0) return NULL;
    if (*i + 1 >= argc) qc_die("option %s needs a value", name);
    *i += 1;
    return argv[*i];
}

static void parse_opts(int argc, char **argv, Opts *o, int *consumed_flags)
{
    (void)consumed_flags;
    for (int i = 0; i < argc; i++) {
        const char *v;
        if ((v = opt_value(argc, argv, &i, "--basis"))) snprintf(o->basis, sizeof o->basis, "%s", v);
        else if ((v = opt_value(argc, argv, &i, "--electrons"))) o->electrons = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--orbitals"))) o->orbitals = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--first"))) o->first = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--maxeval"))) o->maxeval = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--bond"))) o->bond = atof(v);
        else if ((v = opt_value(argc, argv, &i, "--results"))) snprintf(o->results_dir, sizeof o->results_dir, "%s", v);
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) o->verbose = 1;
    }
}

static const char *basis_path(const Opts *o, char *buf, size_t n)
{
    (void)o;
    snprintf(buf, n, "%s/basis/cc-pVQZ.nw", QC_DATA_DIR);
    return buf;
}

static void ensure_ints(Ctx *c, const Opts *o)
{
    if (c->have_ints) return;
    char pth[1024];
    basis_path(o, pth, sizeof pth);
    memset(&c->mol, 0, sizeof c->mol);
    mol_h2_setup(&c->mol, o->bond, 0, 1);
    basis_build(&c->bs, &c->mol, pth, o->basis);
    qc_mkdir_p(CACHE_DIR);
    char cache[512];
    snprintf(cache, sizeof cache, CACHE_DIR "/ao_%s_H2_%.4f.bin", o->basis, o->bond);
    if (ao_integrals_read(&c->ints, cache, c->bs.nao) != 0) {
        double t0 = qc_wall_time();
        ao_integrals_compute(&c->ints, &c->bs, &c->mol);
        printf("  computed AO integrals with libcint in %.2f s\n", qc_wall_time() - t0);
        ao_integrals_write(&c->ints, cache, o->basis);
    } else {
        printf("  loaded AO integrals from %s\n", cache);
    }
    c->have_ints = 1;
}

static void save_rhf(const RHFResult *r, const char *path)
{
    FILE *f = fopen(path, "wb");
    if (!f) return;
    size_t n = (size_t)r->nao * r->nao;
    fwrite(&r->nao, sizeof r->nao, 1, f);
    fwrite(&r->nocc, sizeof r->nocc, 1, f);
    fwrite(&r->energy, sizeof r->energy, 1, f);
    fwrite(&r->enuc, sizeof r->enuc, 1, f);
    fwrite(r->C, sizeof(double), n, f);
    fwrite(r->eps, sizeof(double), (size_t)r->nao, f);
    fwrite(r->P, sizeof(double), n, f);
    fclose(f);
}

static int load_rhf(RHFResult *r, const char *path, int nao)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    memset(r, 0, sizeof *r);
    int ok = fread(&r->nao, sizeof r->nao, 1, f) == 1 && r->nao == nao;
    ok = ok && fread(&r->nocc, sizeof r->nocc, 1, f) == 1;
    ok = ok && fread(&r->energy, sizeof r->energy, 1, f) == 1;
    ok = ok && fread(&r->enuc, sizeof r->enuc, 1, f) == 1;
    size_t n = (size_t)nao * nao;
    r->C = qc_calloc(n, sizeof(double));
    r->eps = qc_calloc((size_t)nao, sizeof(double));
    r->P = qc_calloc(n, sizeof(double));
    ok = ok && fread(r->C, sizeof(double), n, f) == n;
    ok = ok && fread(r->eps, sizeof(double), (size_t)nao, f) == (size_t)nao;
    ok = ok && fread(r->P, sizeof(double), n, f) == n;
    fclose(f);
    if (!ok) { rhf_free(r); return -1; }
    return 0;
}

static void ensure_rhf(Ctx *c, const Opts *o)
{
    ensure_ints(c, o);
    if (c->have_rhf) return;
    char cache[512];
    snprintf(cache, sizeof cache, CACHE_DIR "/rhf_%s_H2_%.4f.bin", o->basis, o->bond);
    if (load_rhf(&c->rhf, cache, c->ints.nao) == 0) {
        printf("  loaded RHF orbitals from %s\n", cache);
    } else {
        double t0 = qc_wall_time();
        int conv = rhf_run(&c->rhf, &c->ints, &c->mol, 200, 1e-10, 1e-8, o->verbose);
        printf("  RHF %s in %d iterations (%.2f s)\n", conv == 0 ? "converged" : "NOT converged",
               c->rhf.iterations, qc_wall_time() - t0);
        save_rhf(&c->rhf, cache);
    }
    c->have_rhf = 1;
}

static int default_first(const Opts *o, const Molecule *mol)
{
    if (o->first >= 0) return o->first;
    return (mol->nelec - o->electrons) / 2;
}

/* active space + Hamiltonian with validation; returns 0 on success */
static int build_active(Ctx *c, const Opts *o)
{
    ensure_rhf(c, o);
    int first = default_first(o, &c->mol);
    active_space_build(&c->as, &c->ints, &c->rhf, &c->mol, o->electrons, o->orbitals, first);
    fermion_build_hamiltonian(&c->fh, &c->as);
    return 0;
}

static void print_rule(void) { printf("----------------------------------------------------------------\n"); }

/* ---------------- commands ---------------- */

int cas_system_load(CasSystem *cs, int nel, int norb, double bond)
{
    Opts o; opts_default(&o);
    o.electrons = nel;
    o.orbitals = norb;
    o.bond = bond;
    Ctx c; memset(&c, 0, sizeof c);
    build_active(&c, &o);
    cs->mol = c.mol;
    cs->rhf = c.rhf;
    cs->as = c.as;
    cs->fh = c.fh;
    ao_integrals_free(&c.ints);
    basis_free(&c.bs);
    return 0;
}

void cas_system_free(CasSystem *cs)
{
    fermion_free(&cs->fh);
    active_space_free(&cs->as);
    rhf_free(&cs->rhf);
}

int cmd_prepare(int argc, char **argv)
{
    Opts o; opts_default(&o); parse_opts(argc, argv, &o, NULL);
    qc_print_banner("prepare");
    Ctx c; memset(&c, 0, sizeof c);
    ensure_ints(&c, &o);
    double maxdev = 0.0;
    for (int i = 0; i < c.ints.nao; i++) {
        double d = fabs(c.ints.S[(size_t)i * c.ints.nao + i] - 1.0);
        if (d > maxdev) maxdev = d;
    }
    printf("  basis            : %s (data/basis/cc-pVQZ.nw)\n", o.basis);
    printf("  nuclear repulsion: %.12f Eh\n", c.mol.enuc);
    printf("  AO functions     : %d (spherical)\n", c.bs.nao);
    printf("  max |S_ii - 1|   : %.3e\n", maxdev);
    double trs = 0, trh = 0;
    for (int i = 0; i < c.ints.nao; i++) { trs += c.ints.S[(size_t)i * c.ints.nao + i]; trh += c.ints.H[(size_t)i * c.ints.nao + i]; }
    printf("  Tr(S) = %.12f ; Tr(H) = %.12f\n", trs, trh);
    if (maxdev > 1e-8) { printf("  FAILED: basis not normalised\n"); return 1; }
    return 0;
}

int cmd_rhf(int argc, char **argv)
{
    Opts o; opts_default(&o); parse_opts(argc, argv, &o, NULL);
    qc_print_banner("RHF (full basis)");
    Ctx c; memset(&c, 0, sizeof c);
    o.verbose = 1;
    ensure_rhf(&c, &o);
    printf("  basis            : %s   (full AO basis, %d functions)\n", o.basis, c.ints.nao);
    printf("  RHF total energy : %.12f Eh\n", c.rhf.energy);
    printf("  nuclear repulsion: %.12f Eh\n", c.mol.enuc);
    printf("  electronic energy: %.12f Eh\n", c.rhf.energy - c.mol.enuc);
    printf("  orbital energies : ");
    for (int i = 0; i < 6 && i < c.ints.nao; i++) printf("%.6f ", c.rhf.eps[i]);
    printf("...\n");
    cJSON *j = cJSON_CreateObject();
    cJSON_AddStringToObject(j, "scope", "full-basis RHF (not active space)");
    cJSON_AddStringToObject(j, "basis", o.basis);
    cJSON_AddNumberToObject(j, "nao", c.ints.nao);
    cJSON_AddNumberToObject(j, "bond_angstrom", o.bond);
    cJSON_AddNumberToObject(j, "rhf_energy_Eh", c.rhf.energy);
    cJSON_AddNumberToObject(j, "nuclear_repulsion_Eh", c.mol.enuc);
    cJSON_AddNumberToObject(j, "iterations", c.rhf.iterations);
    char *meta = results_metadata_json();
    cJSON_AddItemToObject(j, "metadata", cJSON_Parse(meta));
    free(meta);
    char *s = cJSON_Print(j);
    char path[600];
    snprintf(path, sizeof path, "%s/rhf_%s.json", o.results_dir, o.basis);
    results_write_json(path, s);
    free(s);
    cJSON_Delete(j);
    return 0;
}

int cmd_active_space(int argc, char **argv)
{
    Opts o; opts_default(&o); parse_opts(argc, argv, &o, NULL);
    if (2 * o.orbitals > QC_MAX_QUBIT_SV)
        qc_die("CAS(%d,%d) needs %d qubits; the statevector simulator and Pauli cross-checks support at most %d qubits (%d spatial orbitals)",
               o.electrons, o.orbitals, 2 * o.orbitals, QC_MAX_QUBIT_SV, QC_MAX_QUBIT_SV / 2);
    qc_print_banner("active space CAS(2,6)");
    Ctx c; memset(&c, 0, sizeof c);
    build_active(&c, &o);
    printf("  full basis       : %s, %d AO functions (RHF run in the full basis)\n", o.basis, c.ints.nao);
    printf("  RHF energy       : %.12f Eh\n", c.rhf.energy);
    printf("  active space     : CAS(%d,%d) = %d electrons in %d spatial orbitals (MO %d..%d)\n",
           o.electrons, o.orbitals, o.electrons, o.orbitals, c.as.first, c.as.first + o.orbitals - 1);
    printf("  inactive core    : %d doubly occupied MOs folded into E_core\n", c.as.ncore);
    printf("  qubits           : %d (spin orbitals, Jordan-Wigner)\n", c.fh.nspin);
    print_rule();
    double comm, herm, hfe;
    int ok = fermion_validate(&c.fh, &c.as, &c.rhf, &comm, &herm, &hfe, stdout) == 0;
    print_rule();
    printf("  validation       : %s\n", ok ? "PASSED" : "FAILED");

    /* FCI reference (determinant) and Pauli-restricted cross-check */
    FCIResult fci;
    double t0 = qc_wall_time();
    fci_solve(&fci, &c.fh, &c.as, o.electrons);
    double tf = qc_wall_time() - t0;
    double pg = fci_pauli_restricted_ground(&c.fh, &fci);
    printf("  CAS-FCI energy   : %.12f Eh   (%d determinants, %.2f s)\n", fci.ground, fci.ndet, tf);
    printf("  Pauli-restricted : %.12f Eh   |diff| = %.3e\n", pg, fabs(pg - fci.ground));
    printf("  Pauli terms      : %zu\n", c.fh.ham.n);
    printf("  scope            : active-space CAS(%d,%d) energy; NOT the full cc-pVQZ energy\n", o.electrons, o.orbitals);

    /* write results */
    cJSON *j = cJSON_CreateObject();
    cJSON_AddStringToObject(j, "scope", "CAS active space (not full basis)");
    cJSON_AddStringToObject(j, "basis", o.basis);
    cJSON_AddNumberToObject(j, "bond_angstrom", o.bond);
    cJSON_AddNumberToObject(j, "active_electrons", o.electrons);
    cJSON_AddNumberToObject(j, "active_orbitals", o.orbitals);
    cJSON_AddNumberToObject(j, "active_spin_orbitals", c.fh.nspin);
    cJSON_AddNumberToObject(j, "first_mo", c.as.first);
    cJSON_AddNumberToObject(j, "inactive_doubly_occupied", c.as.ncore);
    cJSON_AddNumberToObject(j, "qubits", c.fh.nspin);
    cJSON_AddNumberToObject(j, "pauli_terms", (double)c.fh.ham.n);
    cJSON_AddNumberToObject(j, "rhf_full_basis_Eh", c.rhf.energy);
    cJSON_AddNumberToObject(j, "active_space_hf_Eh", active_hf_energy(&c.as));
    cJSON_AddNumberToObject(j, "fci_cas_Eh", fci.ground);
    cJSON_AddNumberToObject(j, "fci_pauli_restricted_Eh", pg);
    cJSON_AddNumberToObject(j, "validation_passed", ok);
    cJSON_AddNumberToObject(j, "validation_comm_norm", comm);
    cJSON_AddNumberToObject(j, "validation_herm_err", herm);
    cJSON_AddNumberToObject(j, "validation_hf_energy_err", hfe);
    char *meta = results_metadata_json();
    cJSON_AddItemToObject(j, "metadata", cJSON_Parse(meta));
    free(meta);
    char *s = cJSON_Print(j);
    char path[600];
    snprintf(path, sizeof path, "%s/active_CAS%d_%d.json", o.results_dir, o.electrons, o.orbitals);
    results_write_json(path, s);
    free(s);
    cJSON_Delete(j);

    /* Pauli list for reference */
    char lpath[600];
    snprintf(lpath, sizeof lpath, "%s/hamiltonian_CAS%d_%d.txt", o.results_dir, o.electrons, o.orbitals);
    qc_mkdir_p(o.results_dir);
    FILE *f = fopen(lpath, "w");
    if (f) {
        fprintf(f, "# JW Pauli Hamiltonian, qubit 0 leftmost, canonical Y convention\n");
        fprintf(f, "# constant (ecore) term included as identity\n");
        char lab[64];
        for (size_t t = 0; t < c.fh.ham.n; t++) {
            pauli_label(c.fh.ham.t[t].x, c.fh.ham.t[t].z, c.fh.nspin, lab);
            fprintf(f, "%+.15e %+.3e %s\n", creal(c.fh.ham.t[t].c), cimag(c.fh.ham.t[t].c), lab);
        }
        fclose(f);
        printf("  wrote %s\n", lpath);
    }
    fci_free(&fci);
    return ok ? 0 : 2;
}

int cmd_vqe(int argc, char **argv)
{
    Opts o; opts_default(&o);
    const char *mode = "local";
    int submit = 0, yes = 0;
    int nargs = 0;
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) { mode = argv[++i]; continue; }
        if (strcmp(argv[i], "--submit") == 0) { submit = 1; continue; }
        if (strcmp(argv[i], "--yes") == 0) { yes = 1; continue; }
        nargs++;
    }
    (void)nargs;
    parse_opts(argc, argv, &o, NULL);
    if (strcmp(mode, "hardware") == 0) {
        /* hardware path lives in hw/hw_cli.c and needs the same parsed options */
        extern int hw_run_vqe_hardware(int argc, char **argv, int submit, int yes);
        return hw_run_vqe_hardware(argc, argv, submit, yes);
    }
    if (strcmp(mode, "local") != 0) qc_die("unknown --mode '%s' (local|hardware)", mode);

    qc_print_banner("VQE (local statevector, UCCSD + COBYLA)");
    Ctx c; memset(&c, 0, sizeof c);
    build_active(&c, &o);
    printf("  CAS(%d,%d), %d qubits, Pauli terms %zu\n", o.electrons, o.orbitals, c.fh.nspin, c.fh.ham.n);

    UCCSD u;
    uccsd_build(&u, c.fh.nspin, o.electrons);
    printf("  UCCSD: %d spin-conserving excitations (singles+doubles)\n", u.nexc);

    VQEOptions vo = {o.maxeval, 0.1, 1e-9, o.verbose};
    VQEResult vr;
    vqe_optimize_local(&u, &c.fh.ham, &vo, &vr);

    /* references */
    FCIResult fci;
    fci_solve(&fci, &c.fh, &c.as, o.electrons);
    double eref_rhf = c.rhf.energy;

    printf("  COBYLA: %s after %d evaluations (%.2f s)\n", vr.converged ? "converged" : "stopped", vr.nfev, vr.seconds);
    print_rule();
    printf("  RHF (full basis)       : %.12f Eh\n", eref_rhf);
    printf("  CAS(%d,%d) FCI reference : %.12f Eh\n", o.electrons, o.orbitals, fci.ground);
    printf("  CAS(%d,%d) VQE (local)   : %.12f Eh\n", o.electrons, o.orbitals, vr.energy);
    printf("  VQE - FCI              : %+.3e Eh  (variational: must be >= 0)\n", vr.energy - fci.ground);
    print_rule();

    /* logical circuit */
    Circuit circ;
    uccsd_circuit(&u, vr.theta, &circ);
    char qpath[600], jpath[600];
    snprintf(qpath, sizeof qpath, "%s/vqe_local_CAS%d_%d.qasm", o.results_dir, o.electrons, o.orbitals);
    qc_mkdir_p(o.results_dir);
    char *qs = qasm3_string(&circ, circ.nq, 0);
    results_write_json(qpath, qs);
    free(qs);
    printf("  logical circuit: %d gates, depth %d, %d two-qubit gates\n",
           circ.ng, circ_depth(&circ), circ_count_2q(&circ));

    cJSON *j = cJSON_CreateObject();
    cJSON_AddStringToObject(j, "scope", "CAS active space VQE (local statevector)");
    cJSON_AddStringToObject(j, "basis", o.basis);
    cJSON_AddNumberToObject(j, "active_electrons", o.electrons);
    cJSON_AddNumberToObject(j, "active_orbitals", o.orbitals);
    cJSON_AddNumberToObject(j, "qubits", c.fh.nspin);
    cJSON_AddNumberToObject(j, "rhf_full_basis_Eh", eref_rhf);
    cJSON_AddNumberToObject(j, "fci_cas_Eh", fci.ground);
    cJSON_AddNumberToObject(j, "vqe_local_Eh", vr.energy);
    cJSON_AddNumberToObject(j, "vqe_minus_fci_Eh", vr.energy - fci.ground);
    cJSON_AddNumberToObject(j, "nfev", vr.nfev);
    cJSON_AddNumberToObject(j, "converged", vr.converged);
    cJSON_AddNumberToObject(j, "nparam", vr.nparam);
    cJSON *th = cJSON_CreateDoubleArray(vr.theta, vr.nparam);
    cJSON_AddItemToObject(j, "theta", th);
    cJSON_AddNumberToObject(j, "logical_gates", circ.ng);
    cJSON_AddNumberToObject(j, "logical_depth", circ_depth(&circ));
    cJSON_AddNumberToObject(j, "logical_2q", circ_count_2q(&circ));
    char *meta = results_metadata_json();
    cJSON_AddItemToObject(j, "metadata", cJSON_Parse(meta));
    free(meta);
    char *s = cJSON_Print(j);
    snprintf(jpath, sizeof jpath, "%s/vqe_local_CAS%d_%d.json", o.results_dir, o.electrons, o.orbitals);
    results_write_json(jpath, s);
    free(s);
    cJSON_Delete(j);
    circ_free(&circ);
    fci_free(&fci);
    vqe_result_free(&vr);
    uccsd_free(&u);
    (void)submit; (void)yes;
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr,
                "usage: qchem-vqe <command> [options]\n"
                "  prepare --basis cc-pVQZ\n"
                "  rhf [--basis cc-pVQZ]\n"
                "  active-space --electrons 2 --orbitals 6 [--first N]\n"
                "  vqe --mode local [--maxeval N] [--electrons 2 --orbitals 6]\n"
                "  hardware list | inspect [--backend NAME] | preflight [--backend NAME]\n"
                "  vqe --mode hardware [--submit [--yes]]\n"
                "  job status JOB_ID | job results JOB_ID\n");
        return 1;
    }
    const char *cmd = argv[1];
    int rc;
    if (strcmp(cmd, "prepare") == 0) rc = cmd_prepare(argc - 2, argv + 2);
    else if (strcmp(cmd, "rhf") == 0) rc = cmd_rhf(argc - 2, argv + 2);
    else if (strcmp(cmd, "active-space") == 0) rc = cmd_active_space(argc - 2, argv + 2);
    else if (strcmp(cmd, "vqe") == 0) rc = cmd_vqe(argc - 2, argv + 2);
    else if (strcmp(cmd, "hardware") == 0) rc = cmd_hardware(argc - 2, argv + 2);
    else if (strcmp(cmd, "job") == 0) rc = cmd_job(argc - 2, argv + 2);
    else { fprintf(stderr, "unknown command '%s'\n", cmd); rc = 1; }
    return rc;
}
