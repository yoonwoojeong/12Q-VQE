/*
 * qchem-vqe: shared interfaces.
 * All energies are in Hartree, lengths in Bohr unless stated otherwise.
 */
#ifndef QCHEM_H
#define QCHEM_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <complex.h>

#define QC_MAX_ATOM 4
#define QC_MAX_SHELL_PRIM 64
#define QC_MAX_QUBIT_SV 12          /* statevector simulator limit */
#define QC_BOHR_PER_ANG 1.8897261246257702  /* CODATA 2018 */
#define QC_VERSION "0.1.0"
#define QC_PI 3.14159265358979323846

/* ---------- common utilities (src/common/util.c) ---------- */
__attribute__((noreturn)) void qc_die(const char *fmt, ...);
void qc_warn(const char *fmt, ...);
void *qc_calloc(size_t n, size_t sz);
void *qc_realloc(void *p, size_t n);
char *qc_strdup(const char *s);
int qc_mkdir_p(const char *path);
char *qc_read_file(const char *path, size_t *len_out);   /* caller frees */
int qc_write_file(const char *path, const char *data, size_t len);
double qc_wall_time(void);
int qc_env_get(const char *name, char *buf, size_t buflen); /* 1 if set and non-empty */
void qc_print_banner(const char *what);

/* ---------- molecule & basis (src/basis/basis.c) ---------- */
typedef struct {
    int natom;
    int Z[QC_MAX_ATOM];
    double xyz[QC_MAX_ATOM][3];   /* Bohr */
    int charge;
    int mult;
    int nelec;
    double enuc;
} Molecule;

typedef struct {
    int nprim;
    double exp[QC_MAX_SHELL_PRIM];
    double coef[QC_MAX_SHELL_PRIM];   /* contraction coefficients (one spherical contraction) */
    int l;
} Contraction;

typedef struct {
    /* libcint data layout */
    int *atm;      /* 6 * natm */
    int *bas;      /* 8 * nbas */
    double *env;
    int natm;
    int nbas;
    int nenv;
    /* AO bookkeeping */
    int nao;
    int *ao_loc;   /* nbas + 1 entries, offsets of each shell in AO index */
    int *shell_atom;
    int *shell_l;
    char name[64];
} BasisSet;

void mol_h2_setup(Molecule *m, double bond_angstrom, int charge, int mult);
int basis_build(BasisSet *bs, const Molecule *mol, const char *nwchem_path, const char *name);
void basis_free(BasisSet *bs);

/* ---------- AO integrals (src/integrals/integrals.c) ---------- */
typedef struct {
    int nao;
    double *S;    /* nao*nao */
    double *T;
    double *V;
    double *H;    /* T+V */
    double *eri;  /* nao^4, chemist (ij|kl) */
} AOInts;

int ao_integrals_compute(AOInts *ints, const BasisSet *bs, const Molecule *mol);
void ao_integrals_free(AOInts *ints);
int ao_integrals_write(const AOInts *ints, const char *path, const char *tag);
int ao_integrals_read(AOInts *ints, const char *path, int nao_expected);

/* ---------- RHF (src/scf/rhf.c) ---------- */
typedef struct {
    int nao;
    int nocc;          /* doubly occupied orbitals */
    double energy;     /* total, including nuclear repulsion */
    double enuc;
    double *C;         /* nao*nao, columns are MOs (row-major C[mu*nao + p]) */
    double *eps;       /* nao */
    double *P;         /* density, nao*nao */
    int iterations;
    double delta_e;
    double rms_density;
} RHFResult;

int rhf_run(RHFResult *r, const AOInts *ints, const Molecule *mol, int maxiter, double econv, double dconv, int verbose);
void rhf_free(RHFResult *r);

/* ---------- active space (src/active/active.c) ---------- */
typedef struct {
    int nact;          /* spatial active orbitals */
    int nelec_act;     /* active electrons */
    int ncore;         /* doubly occupied inactive orbitals folded into core */
    int first;         /* first MO index of the window */
    double ecore;      /* nuclear repulsion + inactive contribution */
    double *h;         /* nact*nact effective one-body */
    double *g;         /* nact^4 chemist (pq|rs) */
    double *C_act;     /* nao*nact MO coefficients of the active orbitals */
    double *eps_act;   /* nact */
    int nao;
} ActiveSpace;

int active_space_build(ActiveSpace *as, const AOInts *ints, const RHFResult *rhf,
                       const Molecule *mol, int nelec_act, int norb, int first);
double active_hf_energy(const ActiveSpace *as);     /* energy of the active-space HF determinant incl. ecore */
void active_space_free(ActiveSpace *as);

/* ---------- fermion & Pauli algebra (src/fermion, one file per module) ---------- */
typedef struct {
    uint32_t x;        /* X part */
    uint32_t z;        /* Z part */
    double complex c;  /* coefficient of canonical Pauli  i^{pop(x&z)} X^x Z^z */
} PauliTerm;

typedef struct {
    int nq;
    size_t n;
    size_t cap;
    PauliTerm *t;
} PauliSum;

void pauli_sum_init(PauliSum *ps, int nq);
void pauli_sum_free(PauliSum *ps);
void pauli_sum_add(PauliSum *ps, uint32_t x, uint32_t z, double complex c);
void pauli_sum_merge(PauliSum *ps, double tol);           /* combine duplicates, drop |c|<tol */
void pauli_sum_scale(PauliSum *ps, double complex s);
void pauli_sum_append(PauliSum *dst, const PauliSum *src, double complex s);
PauliSum pauli_mul(const PauliSum *a, const PauliSum *b, double tol);
int pauli_label(uint32_t x, uint32_t z, int nq, char *out);   /* out has nq+1 chars, qubit 0 first */
int pauli_is_hermitian(const PauliSum *ps, double tol, double *max_imag);

typedef struct {
    int nspin;               /* 2*norb */
    PauliSum ham;            /* JW Pauli Hamiltonian incl. constant */
    double ecore;
    double *h_spin;          /* nspin^2 */
    double *g_spin;          /* nspin^4, physicist <pq|rs> */
} FermionHam;

int fermion_build_hamiltonian(FermionHam *fh, const ActiveSpace *as);
void fermion_free(FermionHam *fh);
PauliSum jw_number_operator(int nspin);
void fermion_ladder(int j, int dagger, PauliSum *out);   /* JW a_j or a+_j */
int fermion_validate(const FermionHam *fh, const ActiveSpace *as, const RHFResult *rhf,
                     double *comm_norm, double *herm_err, double *hf_energy_err, FILE *log);

/* ---------- statevector simulator (src/sim/statevec.c) ---------- */
typedef struct {
    int nq;
    size_t dim;
    double complex *a;
} StateVec;

int sv_init(StateVec *sv, int nq);
void sv_free(StateVec *sv);
void sv_set_basis(StateVec *sv, size_t index);
void sv_copy(StateVec *dst, const StateVec *src);
void sv_apply_pauli_rotation(StateVec *sv, uint32_t x, uint32_t z, double phi); /* exp(i phi P) */
void sv_apply_1q(StateVec *sv, int q, const double complex u[2][2]);
void sv_apply_cx(StateVec *sv, int c, int t);
void sv_apply_cz(StateVec *sv, int c, int t);
void sv_apply_sum(const StateVec *in, const PauliSum *ps, StateVec *out);
double complex sv_expect_pauli(const StateVec *sv, uint32_t x, uint32_t z);
double sv_expect_sum(const StateVec *sv, const PauliSum *ps);
double sv_norm(const StateVec *sv);

/* ---------- logical / native circuits and OpenQASM 3 (src/qasm/circuit.c) ---------- */
typedef enum {
    G_H, G_SDG, G_S, G_X, G_SX, G_RZ, G_CX, G_CZ, G_MEASURE_NONE
} GateKind;

typedef struct {
    GateKind kind;
    int q0, q1;
    double theta;
} Gate;

typedef struct {
    int nq;          /* register width used in the emitted QASM (physical or logical) */
    int ng, cap;
    Gate *g;
} Circuit;

void circ_init(Circuit *c, int nq);
void circ_free(Circuit *c);
void circ_add(Circuit *c, GateKind k, int q0, int q1, double theta);
void circ_append_pauli_rotation(Circuit *c, uint32_t x, uint32_t z, double phi, int nq);
void circ_to_native(const Circuit *in, Circuit *out);     /* rz, sx, x, cz only */
int circ_depth(const Circuit *c);
int circ_count_2q(const Circuit *c);
int circ_count_kind(const Circuit *c, GateKind k);
void circ_simulate(const Circuit *c, StateVec *sv);
int qasm3_write(const Circuit *c, FILE *f, const char *const *phys_names, int native);
char *qasm3_string(const Circuit *c, int nq_declared, int native);   /* caller frees */

/* ---------- UCCSD ansatz (src/vqe/uccsd.c) ---------- */
typedef struct {
    uint32_t x, z;
    double d;          /* generator coefficient: G = sum_s i*d_s*P_s */
} GenTerm;

typedef struct {
    int kind;          /* 1 single, 2 double */
    int occ[2], vir[2];
    int nterm;
    GenTerm *terms;
} Excitation;

typedef struct {
    int nq;
    int nelec;
    int nexc;
    Excitation *exc;
    uint32_t hf_index;   /* basis index of the HF determinant */
    uint32_t hf_mask;
} UCCSD;

int uccsd_build(UCCSD *u, int nspin, int nelec);
void uccsd_free(UCCSD *u);
void uccsd_prepare(const UCCSD *u, const double *theta, StateVec *sv);   /* from HF */
void uccsd_circuit(const UCCSD *u, const double *theta, Circuit *c);       /* logical IR */

/* ---------- VQE (src/vqe/vqe.c) ---------- */
typedef struct {
    int maxeval;
    double init_step;
    double xtol_rel;
    int verbose;
} VQEOptions;

typedef struct {
    int nparam;
    double *theta;
    double energy;
    int nfev;
    int converged;
    double seconds;
} VQEResult;

int vqe_optimize_local(const UCCSD *u, const PauliSum *ham, const VQEOptions *opt, VQEResult *r);
double vqe_energy(const UCCSD *u, const PauliSum *ham, const double *theta, StateVec *scratch);
void vqe_result_free(VQEResult *r);

/* ---------- FCI (src/fci/fci.c) ---------- */
typedef struct {
    int ndet;
    double *eigvals;      /* all eigenvalues, ascending, including ecore */
    double ground;
    double *ground_vec;   /* ndet */
    uint32_t *dets;       /* occupation bitmasks over spin orbitals */
} FCIResult;

int fci_solve(FCIResult *r, const FermionHam *fh, const ActiveSpace *as, int nelec);
double fci_pauli_restricted_ground(const FermionHam *fh, const FCIResult *fci);
void fci_free(FCIResult *r);

/* ---------- results / metadata (src/common/results.c) ---------- */
int results_write_json(const char *path, const char *json_text);
char *results_metadata_json(void);   /* caller frees */

/* ---------- IBM Quantum REST (src/hw, one file per module) ---------- */
typedef struct {
    char api_key[256];
    char crn[512];
    char base_url[256];
    char api_version[32];
    char token[4096];
    double token_expiry;      /* seconds since epoch */
} IBMSession;

typedef struct {
    char name[128];
    int n_qubits;
    int is_simulator;
    char status[32];
    int queue_length;
} IBMBackendInfo;

#define IBM_MAX_QUBIT 256
#define IBM_MAX_EDGE 1024

typedef struct {
    int nqubit;
    int nedge;                       /* undirected couplers */
    int edge_a[IBM_MAX_EDGE], edge_b[IBM_MAX_EDGE];
    double edge_err[IBM_MAX_EDGE];   /* CZ error, NaN if unknown */
    double edge_len_ns[IBM_MAX_EDGE];
    double readout_err[IBM_MAX_QUBIT];
    double sx_err[IBM_MAX_QUBIT];
    double x_err[IBM_MAX_QUBIT];
    double t1_us[IBM_MAX_QUBIT];
    double t2_us[IBM_MAX_QUBIT];
    int has_cz, has_ecr;
    char backend_name[128];
    char backend_version[32];
} IBMCalibration;

typedef struct {
    int nl;                          /* logical qubits */
    int phys[QC_MAX_QUBIT_SV];       /* initial logical -> physical */
    int final_pos[QC_MAX_QUBIT_SV];  /* after routing */
    int subset[IBM_MAX_QUBIT];       /* physical qubits in the layout */
    int nsubset;
    double score;                    /* summed error of the chosen subgraph (lower is better) */
    int nswaps;
    int ok;
} IBMLayout;

int ibm_session_env(IBMSession *s);                /* stored file, then env overrides; never prints secrets */
int cred_default_path(char *buf, size_t n);
int cred_load(const char *path, IBMSession *s);
int cred_store(const char *path, const char *api_key, const char *crn, const char *url, const char *version);
void cred_mask(const char *secret, char *out, size_t n);
int ibm_crn_region(const char *crn, char *out, size_t n);
const char *ibm_url_for_crn(const char *crn);
void ibm_report_http_error(const IBMSession *s, const char *method, const char *path, long code, const char *resp);
int cmd_auth(int argc, char **argv);
int ibm_session_authenticate(IBMSession *s);       /* IAM API-key exchange */
int ibm_http(IBMSession *s, const char *method, const char *path, const char *body, char **resp, long *status);
int ibm_list_backends(IBMSession *s, IBMBackendInfo **out, int *nout);
int ibm_get_backend_raw(IBMSession *s, const char *name, char **cfg_json, char **props_json, char **status_json);
int ibm_calibration_parse(const char *cfg_json, const char *props_json, IBMCalibration *cal);
int ibm_job_submit(IBMSession *s, const char *backend, const char *qasm, const char *const *labels, int nlabels,
                   double precision, char **resp);
int ibm_job_get(IBMSession *s, const char *job_id, char **json);
int ibm_job_results(IBMSession *s, const char *job_id, char **json);
int base64_decode_f64(const char *b64, double **out, size_t *n);

int layout_select(const IBMCalibration *cal, int nlogical, IBMLayout *out);
int layout_compile(const Circuit *logical, const IBMCalibration *cal, IBMLayout *lay, Circuit *native_phys);

/* ---------- CAS system loader (src/main.c) ---------- */
typedef struct {
    Molecule mol;
    RHFResult rhf;
    ActiveSpace as;
    FermionHam fh;
} CasSystem;
int cas_system_load(CasSystem *cs, int nel, int norb, double bond);
void cas_system_free(CasSystem *cs);

/* ---------- CLI (src/main.c) ---------- */
int cmd_prepare(int argc, char **argv);
int cmd_rhf(int argc, char **argv);
int cmd_active_space(int argc, char **argv);
int cmd_vqe(int argc, char **argv);
int cmd_hardware(int argc, char **argv);
int cmd_job(int argc, char **argv);

#endif
