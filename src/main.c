/*
 * qchem-vqe command-line interface.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cjson/cJSON.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#ifndef QC_DATA_DIR
#define QC_DATA_DIR "data"
#endif

typedef struct {
    char mol[64];
    char basis[64];
    char xyz_path[512];
    char fcidump_in[512];
    char fcidump_out[512];
    char results_dir[512];
    char cache_dir[512];
    int charge;
    int mult;
    int electrons;      /* active electrons */
    int orbitals;       /* active spatial orbitals */
    int first;          /* first active MO, -1 = default */
    int maxeval;
    int verbose;
    int use_cobyla;
    int use_adapt;
    int adapt_max;
    double adapt_tol;
    double bond;        /* Angstrom */
    double rmin, rmax;  /* PES scan */
    int npts;
    int shots;          /* SQD shots */
    double noise;       /* SQD per-qubit bit-flip noise */
    int sc_iters;       /* SQD self-consistent recovery iterations */
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
    snprintf(o->mol, sizeof o->mol, "H2");
    snprintf(o->basis, sizeof o->basis, "cc-pVQZ");
    qc_results_dir(o->results_dir, sizeof o->results_dir);
    qc_cache_dir(o->cache_dir, sizeof o->cache_dir);
    o->charge = 0;
    o->mult = 1;
    o->electrons = 2;
    o->orbitals = 6;
    o->first = -1;
    o->maxeval = 4000;
    o->use_cobyla = 0;   /* default: L-BFGS with analytic adjoint gradients */
    o->use_adapt = 0;
    o->adapt_max = 35;
    o->adapt_tol = 1e-4;
    o->bond = 0.735;
    o->rmin = 0.4;
    o->rmax = 2.8;
    o->npts = 13;
    o->shots = 4096;
    o->noise = 0.02;
    o->sc_iters = 4;
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
        if ((v = opt_value(argc, argv, &i, "--mol"))) snprintf(o->mol, sizeof o->mol, "%s", v);
        else if ((v = opt_value(argc, argv, &i, "--basis"))) snprintf(o->basis, sizeof o->basis, "%s", v);
        else if ((v = opt_value(argc, argv, &i, "--electrons"))) o->electrons = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--orbitals"))) o->orbitals = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--first"))) o->first = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--maxeval"))) o->maxeval = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--bond"))) o->bond = atof(v);
        else if ((v = opt_value(argc, argv, &i, "--rmin"))) o->rmin = atof(v);
        else if ((v = opt_value(argc, argv, &i, "--rmax"))) o->rmax = atof(v);
        else if ((v = opt_value(argc, argv, &i, "--npts"))) o->npts = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--results"))) snprintf(o->results_dir, sizeof o->results_dir, "%s", v);
        else if ((v = opt_value(argc, argv, &i, "--opt"))) {
            if (strcmp(v, "cobyla") == 0) o->use_cobyla = 1;
            else if (strcmp(v, "lbfgs") == 0) o->use_cobyla = 0;
            else qc_die("--opt must be 'lbfgs' or 'cobyla'");
        }
        else if ((v = opt_value(argc, argv, &i, "--ansatz"))) {
            if (strcmp(v, "adapt") == 0) o->use_adapt = 1;
            else if (strcmp(v, "uccsd") == 0 || strcmp(v, "givens") == 0) o->use_adapt = 0;
        }
        else if ((v = opt_value(argc, argv, &i, "--adapt-max"))) o->adapt_max = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--adapt-tol"))) o->adapt_tol = atof(v);
        else if ((v = opt_value(argc, argv, &i, "--xyz"))) snprintf(o->xyz_path, sizeof o->xyz_path, "%s", v);
        else if ((v = opt_value(argc, argv, &i, "--charge"))) o->charge = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--mult"))) o->mult = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--fcidump"))) snprintf(o->fcidump_in, sizeof o->fcidump_in, "%s", v);
        else if ((v = opt_value(argc, argv, &i, "--fcidump-out"))) snprintf(o->fcidump_out, sizeof o->fcidump_out, "%s", v);
        else if ((v = opt_value(argc, argv, &i, "--shots"))) o->shots = atoi(v);
        else if ((v = opt_value(argc, argv, &i, "--noise"))) o->noise = atof(v);
        else if ((v = opt_value(argc, argv, &i, "--sc-iters"))) o->sc_iters = atoi(v);
        else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) o->verbose = 1;
    }
}

static const char *basis_path(const Opts *o, char *buf, size_t n)
{
    char ddir[512];
    qc_data_dir(ddir, sizeof ddir);
    if (strcmp(o->basis, "sto-3g") == 0 || strcmp(o->basis, "STO-3G") == 0)
        snprintf(buf, n, "%s/basis/sto-3g.nw", ddir);
    else if (strcmp(o->basis, "cc-pVDZ") == 0 || strcmp(o->basis, "cc-pvdz") == 0)
        snprintf(buf, n, "%s/basis/cc-pVDZ.nw", ddir);
    else if (strcmp(o->basis, "cc-pVQZ") == 0 || strcmp(o->basis, "cc-pvqz") == 0)
        snprintf(buf, n, "%s/basis/cc-pVQZ.nw", ddir);
    else
        snprintf(buf, n, "%s/basis/%s.nw", ddir, o->basis);
    return buf;
}

static void ensure_ints(Ctx *c, const Opts *o)
{
    if (c->have_ints) return;
    char pth[1024];
    basis_path(o, pth, sizeof pth);
    memset(&c->mol, 0, sizeof c->mol);
    if (o->xyz_path[0]) {
        if (mol_load_xyz(&c->mol, o->xyz_path, o->charge, o->mult) != 0)
            qc_die("cannot load XYZ geometry from '%s'", o->xyz_path);
        basis_build(&c->bs, &c->mol, pth, o->basis);
        double t0 = qc_wall_time();
        ao_integrals_compute(&c->ints, &c->bs, &c->mol);
        printf("  computed AO integrals from %s with libcint in %.2f s\n", o->xyz_path, qc_wall_time() - t0);
        c->have_ints = 1;
        return;
    }
    if (mol_setup_by_name(&c->mol, o->mol, o->bond) != 0)
        qc_die("unknown molecule '%s' (supported: H2, HeH+, LiH, BeH2, H4 or --xyz FILE)", o->mol);
    basis_build(&c->bs, &c->mol, pth, o->basis);
    qc_mkdir_p(o->cache_dir);
    char cache[1024];
    snprintf(cache, sizeof cache, "%s/ao_%s_%s_%.4f.bin", o->cache_dir, o->basis, o->mol, o->bond);
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
    if (o->xyz_path[0]) {
        double t0 = qc_wall_time();
        int conv = rhf_run(&c->rhf, &c->ints, &c->mol, 200, 1e-10, 1e-8, o->verbose);
        printf("  RHF %s in %d iterations (%.2f s)\n", conv == 0 ? "converged" : "NOT converged",
               c->rhf.iterations, qc_wall_time() - t0);
        c->have_rhf = 1;
        return;
    }
    char cache[1024];
    snprintf(cache, sizeof cache, "%s/rhf_%s_%s_%.4f.bin", o->cache_dir, o->basis, o->mol, o->bond);
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
static int build_active(Ctx *c, Opts *o)
{
    if (o->fcidump_in[0]) {
        if (fcidump_read(&c->as, o->fcidump_in) != 0)
            qc_die("cannot read FCIDUMP file '%s'", o->fcidump_in);
        o->electrons = c->as.nelec_act;
        o->orbitals = c->as.nact;
        c->rhf.energy = active_hf_energy(&c->as);
        fermion_build_hamiltonian(&c->fh, &c->as);
        if (o->fcidump_out[0]) fcidump_write(&c->as, o->fcidump_out);
        return 0;
    }
    ensure_rhf(c, o);
    int first = default_first(o, &c->mol);
    active_space_build(&c->as, &c->ints, &c->rhf, &c->mol, o->electrons, o->orbitals, first);
    fermion_build_hamiltonian(&c->fh, &c->as);
    if (o->fcidump_out[0]) {
        if (fcidump_write(&c->as, o->fcidump_out) == 0)
            printf("  exported FCIDUMP to %s\n", o->fcidump_out);
    }
    return 0;
}

static void print_rule(void) { printf("----------------------------------------------------------------\n"); }

/* ---------------- commands ---------------- */

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
    printf("  molecule         : %s (%d electrons, R = %.3f A)\n", o.mol, c.mol.nelec, o.bond);
    printf("  basis            : %s\n", o.basis);
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
    printf("  molecule         : %s (R = %.3f A)\n", o.mol, o.bond);
    printf("  basis            : %s   (full AO basis, %d functions)\n", o.basis, c.ints.nao);
    printf("  RHF total energy : %.12f Eh\n", c.rhf.energy);
    printf("  nuclear repulsion: %.12f Eh\n", c.mol.enuc);
    printf("  electronic energy: %.12f Eh\n", c.rhf.energy - c.mol.enuc);
    printf("  orbital energies : ");
    for (int i = 0; i < 6 && i < c.ints.nao; i++) printf("%.6f ", c.rhf.eps[i]);
    printf("...\n");
    cJSON *j = cJSON_CreateObject();
    cJSON_AddStringToObject(j, "scope", "full-basis RHF (not active space)");
    cJSON_AddStringToObject(j, "molecule", o.mol);
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
        qc_die("CAS(%d,%d) needs %d qubits; the statevector simulator supports at most %d qubits (%d spatial orbitals)",
               o.electrons, o.orbitals, 2 * o.orbitals, QC_MAX_QUBIT_SV, QC_MAX_QUBIT_SV / 2);
    char banner[128];
    snprintf(banner, sizeof banner, "active space %s CAS(%d,%d)", o.mol, o.electrons, o.orbitals);
    qc_print_banner(banner);
    Ctx c; memset(&c, 0, sizeof c);
    build_active(&c, &o);
    printf("  molecule         : %s (%d total electrons, R = %.3f A)\n", o.mol, c.mol.nelec, o.bond);
    printf("  full basis       : %s, %d AO functions (RHF run in the full basis)\n", o.basis, c.ints.nao);
    printf("  RHF energy       : %.12f Eh\n", c.rhf.energy);
    printf("  active space     : CAS(%d,%d) = %d electrons in %d spatial orbitals (MO %d..%d)\n",
           o.electrons, o.orbitals, o.electrons, o.orbitals, c.as.first, c.as.first + o.orbitals - 1);
    printf("  inactive core    : %d doubly occupied MOs folded into E_core (%.8f Eh)\n", c.as.ncore, c.as.ecore);
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
    printf("  scope            : active-space CAS(%d,%d) energy; NOT the full %s energy\n", o.electrons, o.orbitals, o.basis);

    /* write results */
    cJSON *j = cJSON_CreateObject();
    cJSON_AddStringToObject(j, "scope", "CAS active space (not full basis)");
    cJSON_AddStringToObject(j, "molecule", o.mol);
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
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) { mode = argv[++i]; continue; }
        if (strcmp(argv[i], "--submit") == 0) { submit = 1; continue; }
        if (strcmp(argv[i], "--yes") == 0) { yes = 1; continue; }
    }
    parse_opts(argc, argv, &o, NULL);
    if (strcmp(mode, "hardware") == 0) {
        extern int hw_run_vqe_hardware(int argc, char **argv, int submit, int yes);
        return hw_run_vqe_hardware(argc, argv, submit, yes);
    }
    if (strcmp(mode, "local") != 0) qc_die("unknown --mode '%s' (local|hardware)", mode);

    char banner[128];
    snprintf(banner, sizeof banner, "VQE (%s, %s)", o.use_adapt ? "ADAPT-VQE + L-BFGS" : "UCCSD",
             o.use_cobyla ? "COBYLA" : "L-BFGS + Adjoint Grad");
    qc_print_banner(banner);
    Ctx c; memset(&c, 0, sizeof c);
    build_active(&c, &o);
    printf("  %s (%s, R=%.3f A), CAS(%d,%d), %d qubits, Pauli terms %zu\n",
           o.mol, o.basis, o.bond, o.electrons, o.orbitals, c.fh.nspin, c.fh.ham.n);

    UCCSD u;
    uccsd_build(&u, c.fh.nspin, o.electrons);
    printf("  operator pool: %d spin-conserving excitations (singles+doubles)\n", u.nexc);

    FCIResult fci;
    fci_solve(&fci, &c.fh, &c.as, o.electrons);
    double eref_rhf = c.rhf.energy;

    Circuit circ;
    double vqe_E = 0.0, secs = 0.0;
    int nfev = 0, conv = 0, nparam = 0;
    double *theta_save = NULL;

    if (o.use_adapt) {
        AdaptResult ar;
        adapt_vqe_optimize(&u, &c.fh.ham, o.adapt_max, o.adapt_tol, o.verbose, &ar, &circ);
        vqe_E = ar.energy;
        secs = ar.seconds;
        nfev = ar.nfev;
        conv = ar.converged;
        nparam = ar.n_ops;
        theta_save = qc_calloc((size_t)(nparam ? nparam : 1), sizeof(double));
        if (nparam) memcpy(theta_save, ar.theta, (size_t)nparam * sizeof(double));
        printf("  ADAPT-VQE: %s with %d operators (%d L-BFGS evals, %.3f s), max|g| = %.2e\n",
               conv ? "converged" : "stopped", ar.n_ops, ar.nfev, ar.seconds, ar.max_grad);
        adapt_result_free(&ar);
    } else {
        VQEOptions vo = {o.maxeval, 0.1, 1e-10, o.verbose, o.use_cobyla};
        VQEResult vr;
        vqe_optimize_local(&u, &c.fh.ham, &vo, &vr);
        vqe_E = vr.energy;
        secs = vr.seconds;
        nfev = vr.nfev;
        conv = vr.converged;
        nparam = vr.nparam;
        theta_save = qc_calloc((size_t)nparam, sizeof(double));
        memcpy(theta_save, vr.theta, (size_t)nparam * sizeof(double));
        uccsd_circuit(&u, vr.theta, &circ);
        printf("  %s: %s after %d evaluations (%.3f s)\n",
               o.use_cobyla ? "COBYLA" : "L-BFGS", vr.converged ? "converged" : "stopped", vr.nfev, vr.seconds);
        vqe_result_free(&vr);
    }

    print_rule();
    printf("  RHF (full basis)       : %.12f Eh\n", eref_rhf);
    printf("  CAS(%d,%d) FCI reference : %.12f Eh\n", o.electrons, o.orbitals, fci.ground);
    printf("  CAS(%d,%d) VQE (local)   : %.12f Eh\n", o.electrons, o.orbitals, vqe_E);
    printf("  VQE - FCI              : %+.3e Eh  (variational: must be >= 0)\n", vqe_E - fci.ground);
    print_rule();

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
    cJSON_AddStringToObject(j, "molecule", o.mol);
    cJSON_AddStringToObject(j, "basis", o.basis);
    cJSON_AddStringToObject(j, "ansatz", o.use_adapt ? "adapt" : "uccsd");
    cJSON_AddStringToObject(j, "optimizer", o.use_cobyla ? "cobyla" : "lbfgs");
    cJSON_AddNumberToObject(j, "active_electrons", o.electrons);
    cJSON_AddNumberToObject(j, "active_orbitals", o.orbitals);
    cJSON_AddNumberToObject(j, "qubits", c.fh.nspin);
    cJSON_AddNumberToObject(j, "rhf_full_basis_Eh", eref_rhf);
    cJSON_AddNumberToObject(j, "fci_cas_Eh", fci.ground);
    cJSON_AddNumberToObject(j, "vqe_local_Eh", vqe_E);
    cJSON_AddNumberToObject(j, "vqe_minus_fci_Eh", vqe_E - fci.ground);
    cJSON_AddNumberToObject(j, "nfev", nfev);
    cJSON_AddNumberToObject(j, "seconds", secs);
    cJSON_AddNumberToObject(j, "converged", conv);
    cJSON_AddNumberToObject(j, "nparam", nparam);
    cJSON *th = cJSON_CreateDoubleArray(theta_save, nparam);
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
    free(theta_save);
    circ_free(&circ);
    fci_free(&fci);
    uccsd_free(&u);
    fermion_free(&c.fh);
    active_space_free(&c.as);
    rhf_free(&c.rhf);
    ao_integrals_free(&c.ints);
    basis_free(&c.bs);
    return 0;
}

int cmd_pes(int argc, char **argv)
{
    Opts o; opts_default(&o);
    parse_opts(argc, argv, &o, NULL);
    char banner[256];
    snprintf(banner, sizeof banner, "Potential Energy Surface Scan: %s (%s, CAS(%d,%d))",
             o.mol, o.basis, o.electrons, o.orbitals);
    qc_print_banner(banner);

    cJSON *root = cJSON_CreateObject();
    cJSON_AddStringToObject(root, "molecule", o.mol);
    cJSON_AddStringToObject(root, "basis", o.basis);
    cJSON_AddNumberToObject(root, "active_electrons", o.electrons);
    cJSON_AddNumberToObject(root, "active_orbitals", o.orbitals);
    cJSON *pts = cJSON_AddArrayToObject(root, "points");

    printf("  %-8s  %-16s  %-16s  %-16s  %-16s  %-10s  %-8s\n",
           "R (A)", "RHF (Eh)", "CAS-FCI (Eh)", "UCCSD-VQE (Eh)", "ADAPT-VQE (Eh)", "|VQE-FCI|", "ADAPT 2Q");
    print_rule();

    for (int i = 0; i < o.npts; i++) {
        double r = (o.npts == 1) ? o.bond : o.rmin + (o.rmax - o.rmin) * (double)i / (double)(o.npts - 1);
        Opts oi = o;
        oi.bond = r;
        oi.verbose = 0;
        Ctx c; memset(&c, 0, sizeof c);

        /* silence per-point cache prints by building directly */
        char pth[1024];
        basis_path(&oi, pth, sizeof pth);
        mol_setup_by_name(&c.mol, oi.mol, r);
        basis_build(&c.bs, &c.mol, pth, oi.basis);
        ao_integrals_compute(&c.ints, &c.bs, &c.mol);
        rhf_run(&c.rhf, &c.ints, &c.mol, 200, 1e-10, 1e-8, 0);
        int first = default_first(&oi, &c.mol);
        active_space_build(&c.as, &c.ints, &c.rhf, &c.mol, oi.electrons, oi.orbitals, first);
        fermion_build_hamiltonian(&c.fh, &c.as);

        FCIResult fci;
        fci_solve(&fci, &c.fh, &c.as, oi.electrons);

        UCCSD u;
        uccsd_build(&u, c.fh.nspin, oi.electrons);
        VQEOptions vo = {1000, 0.1, 1e-10, 0, 0};
        VQEResult vr;
        vqe_optimize_local(&u, &c.fh.ham, &vo, &vr);

        AdaptResult ar;
        Circuit acirc;
        adapt_vqe_optimize(&u, &c.fh.ham, oi.adapt_max, oi.adapt_tol, 0, &ar, &acirc);
        int adapt_2q = circ_count_2q(&acirc);

        printf("  %-8.3f  %16.10f  %16.10f  %16.10f  %16.10f  %10.2e  %4d (%2d ops)\n",
               r, c.rhf.energy, fci.ground, vr.energy, ar.energy, fabs(vr.energy - fci.ground), adapt_2q, ar.n_ops);

        cJSON *pt = cJSON_CreateObject();
        cJSON_AddNumberToObject(pt, "bond_angstrom", r);
        cJSON_AddNumberToObject(pt, "rhf_Eh", c.rhf.energy);
        cJSON_AddNumberToObject(pt, "fci_Eh", fci.ground);
        cJSON_AddNumberToObject(pt, "uccsd_vqe_Eh", vr.energy);
        cJSON_AddNumberToObject(pt, "adapt_vqe_Eh", ar.energy);
        cJSON_AddNumberToObject(pt, "adapt_ops", ar.n_ops);
        cJSON_AddNumberToObject(pt, "adapt_2q_gates", adapt_2q);
        cJSON_AddItemToArray(pts, pt);

        circ_free(&acirc);
        adapt_result_free(&ar);
        vqe_result_free(&vr);
        uccsd_free(&u);
        fci_free(&fci);
        fermion_free(&c.fh);
        active_space_free(&c.as);
        rhf_free(&c.rhf);
        ao_integrals_free(&c.ints);
        basis_free(&c.bs);
    }
    print_rule();
    char jpath[1024];
    snprintf(jpath, sizeof jpath, "%s/pes_%s_%s_CAS%d_%d.json", o.results_dir, o.mol, o.basis, o.electrons, o.orbitals);
    char *s = cJSON_Print(root);
    results_write_json(jpath, s);
    free(s);
    cJSON_Delete(root);
    return 0;
}

int cmd_sqd(int argc, char **argv)
{
    Opts o; opts_default(&o);
    parse_opts(argc, argv, &o, NULL);
    qc_print_banner("Sample-Based Quantum Diagonalization (SQD / QSCI)");
    Ctx c; memset(&c, 0, sizeof c);
    build_active(&c, &o);

    FCIResult fci;
    fci_solve(&fci, &c.fh, &c.as, o.electrons);

    /* Prepare quantum state via ADAPT-VQE (or UCCSD) to sample from */
    UCCSD u;
    uccsd_build(&u, c.fh.nspin, o.electrons);
    AdaptResult ar;
    Circuit acirc;
    adapt_vqe_optimize(&u, &c.fh.ham, o.adapt_max, o.adapt_tol, 0, &ar, &acirc);
    StateVec psi;
    sv_init(&psi, c.fh.nspin);
    circ_simulate(&acirc, &psi);

    SQDResult sqd;
    double t0 = qc_wall_time();
    sqd_run_simulated(&c.fh, &c.as, &psi, o.shots, o.noise, o.sc_iters, 42, &sqd);
    double t_sqd = qc_wall_time() - t0;

    printf("  system                     : %s (%s), CAS(%d,%d), %d qubits\n",
           o.fcidump_in[0] ? o.fcidump_in : o.mol, o.basis, o.electrons, o.orbitals, c.fh.nspin);
    printf("  shots & noise              : %d shots, %.2f%% per-qubit bit-flip error, %d SC recovery iterations\n",
           o.shots, 100.0 * o.noise, o.sc_iters);
    printf("  sector retention           : %d / %d shots in (N_a,N_b) sector (%d recovered by S-CORE)\n",
           sqd.n_valid_sector, sqd.n_shots, sqd.n_recovered);
    printf("  sampled subspace dimension : %d / %d total CAS determinants (%.2f ms)\n",
           sqd.subspace_dim, fci.ndet, 1e3 * t_sqd);
    print_rule();
    printf("  RHF reference              : %.12f Eh\n", c.rhf.energy);
    printf("  CAS-FCI exact reference    : %.12f Eh\n", fci.ground);
    printf("  Noisy Estimator <H>        : %.12f Eh   (error %+.6f Eh)\n",
           sqd.noisy_estimator_energy, sqd.noisy_estimator_energy - fci.ground);
    printf("  SQD recovered energy       : %.12f Eh   (error %+.3e Eh)\n",
           sqd.energy, sqd.energy - fci.ground);
    print_rule();

    cJSON *j = cJSON_CreateObject();
    cJSON_AddStringToObject(j, "scope", "Sample-Based Quantum Diagonalization (SQD) with configuration recovery");
    cJSON_AddStringToObject(j, "molecule", o.mol);
    cJSON_AddStringToObject(j, "basis", o.basis);
    cJSON_AddNumberToObject(j, "active_electrons", o.electrons);
    cJSON_AddNumberToObject(j, "active_orbitals", o.orbitals);
    cJSON_AddNumberToObject(j, "qubits", c.fh.nspin);
    cJSON_AddNumberToObject(j, "shots", o.shots);
    cJSON_AddNumberToObject(j, "bitflip_prob", o.noise);
    cJSON_AddNumberToObject(j, "sc_recovery_iters", o.sc_iters);
    cJSON_AddNumberToObject(j, "valid_sector_shots", sqd.n_valid_sector);
    cJSON_AddNumberToObject(j, "recovered_shots", sqd.n_recovered);
    cJSON_AddNumberToObject(j, "subspace_dim", sqd.subspace_dim);
    cJSON_AddNumberToObject(j, "fci_dim", fci.ndet);
    cJSON_AddNumberToObject(j, "rhf_Eh", c.rhf.energy);
    cJSON_AddNumberToObject(j, "fci_Eh", fci.ground);
    cJSON_AddNumberToObject(j, "noisy_estimator_Eh", sqd.noisy_estimator_energy);
    cJSON_AddNumberToObject(j, "sqd_energy_Eh", sqd.energy);
    cJSON_AddNumberToObject(j, "sqd_minus_fci_Eh", sqd.energy - fci.ground);
    char *s = cJSON_Print(j);
    char jpath[1024];
    snprintf(jpath, sizeof jpath, "%s/sqd_%s_CAS%d_%d.json", o.results_dir, o.mol, o.electrons, o.orbitals);
    results_write_json(jpath, s);
    free(s);
    cJSON_Delete(j);

    sqd_result_free(&sqd);
    sv_free(&psi);
    circ_free(&acirc);
    adapt_result_free(&ar);
    uccsd_free(&u);
    fci_free(&fci);
    fermion_free(&c.fh);
    active_space_free(&c.as);
    if (c.have_rhf) rhf_free(&c.rhf);
    if (c.have_ints) { ao_integrals_free(&c.ints); basis_free(&c.bs); }
    return 0;
}

/* Under `wsl --exec` no shell reports a crash, so say it ourselves (async-signal-safe). */
static void crash_handler(int sig)
{
    static const char msg[] = "\nerror: qchem-vqe crashed (internal error). Please report the command and its output.\n";
    ssize_t w = write(STDERR_FILENO, msg, sizeof msg - 1);
    (void)w;
    _exit(128 + sig);
}

int main(int argc, char **argv)
{
    signal(SIGSEGV, crash_handler);
    signal(SIGBUS, crash_handler);
    signal(SIGFPE, crash_handler);
    signal(SIGABRT, crash_handler);
    if (argc < 2) {
        fprintf(stderr,
                "usage: qchem-vqe <command> [options]\n"
                "  prepare [--mol H2|HeH+|LiH|BeH2|H4 | --xyz FILE] [--basis cc-pVQZ|cc-pVDZ|sto-3g] [--bond R]\n"
                "  rhf [--mol NAME | --xyz FILE] [--basis NAME] [--bond R]\n"
                "  active-space --electrons N --orbitals M [--first K] [--fcidump FILE] [--fcidump-out FILE]\n"
                "  vqe --mode local [--ansatz uccsd|adapt] [--opt lbfgs|cobyla] [--mol NAME | --fcidump FILE]\n"
                "  sqd [--mol NAME | --fcidump FILE] [--shots N] [--noise P] [--sc-iters K]\n"
                "  pes [--mol NAME] [--basis NAME] [--rmin R0 --rmax R1 --npts N]\n"
                "  hardware list | inspect [--backend NAME] | preflight [--backend NAME]\n"
                "  vqe --mode hardware [--submit [--yes]]\n"
                "  job status JOB_ID | job results JOB_ID\n"
                "  auth login | auth status | auth logout     IBM credentials (hidden prompt)\n"
                "  experiment plan|run [CONFIG] | experiment results JOB_ID   controlled QPU experiment\n");
        return 1;
    }
    const char *cmd = argv[1];
    int rc;
    if (strcmp(cmd, "prepare") == 0) rc = cmd_prepare(argc - 2, argv + 2);
    else if (strcmp(cmd, "rhf") == 0) rc = cmd_rhf(argc - 2, argv + 2);
    else if (strcmp(cmd, "active-space") == 0) rc = cmd_active_space(argc - 2, argv + 2);
    else if (strcmp(cmd, "vqe") == 0) rc = cmd_vqe(argc - 2, argv + 2);
    else if (strcmp(cmd, "sqd") == 0) rc = cmd_sqd(argc - 2, argv + 2);
    else if (strcmp(cmd, "pes") == 0) rc = cmd_pes(argc - 2, argv + 2);
    else if (strcmp(cmd, "hardware") == 0) rc = cmd_hardware(argc - 2, argv + 2);
    else if (strcmp(cmd, "job") == 0) rc = cmd_job(argc - 2, argv + 2);
    else if (strcmp(cmd, "auth") == 0) rc = cmd_auth(argc - 2, argv + 2);
    else if (strcmp(cmd, "experiment") == 0) rc = cmd_experiment(argc - 2, argv + 2);
    else { fprintf(stderr, "unknown command '%s'\n", cmd); rc = 1; }
    return rc;
}
