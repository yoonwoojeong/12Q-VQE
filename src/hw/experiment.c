/*
 * Controlled QPU experiment, defined by one small config file and driven by three commands:
 *
 *   qchem-vqe experiment plan    [CONFIG] [--calibration CONF.json PROPS.json]
 *   qchem-vqe experiment run     [CONFIG] [--yes]
 *   qchem-vqe experiment results JOB_ID
 *
 * All arms go into ONE Estimator job (one PUB per arm), so they share the backend, the
 * physical qubits, the calibration window, shots, mitigation and seed:
 *
 *   bare_hf     X gates only                         -> readout-only control, exact E = E_HF
 *   matched_hf  the VQE circuit with every rotation  -> same gates and noise as the treatment,
 *               angle set to 0 (prepares HF exactly)    exact E = E_HF
 *   vqe         the circuit preparing the VQE state  -> treatment, exact E = noiseless circuit energy
 *
 * The matched control measures the noise bias of the treatment circuit; the correlation
 * energy is measured as E(vqe) - E(matched_hf) on the same chip in the same job.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cjson/cJSON.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define EXP_DIR "results/experiments"
#define EXP_DEFAULT_CONFIG "experiments/h2_cas26.ini"

/* optimised UCCSD parameters, cached per experiment name and checked against the config */
static double *exp_theta(const ExpConfig *cfg, const UCCSD *u, const PauliSum *ham, double *e_vqe)
{
    char path[512];
    snprintf(path, sizeof path, EXP_DIR "/%s_vqe.json", cfg->name);
    double *theta = qc_calloc((size_t)u->nexc, sizeof(double));
    size_t len = 0;
    char *txt = qc_read_file(path, &len);
    if (txt) {
        cJSON *j = cJSON_Parse(txt);
        free(txt);
        cJSON *th = cJSON_GetObjectItem(j, "theta");
        cJSON *b = cJSON_GetObjectItem(j, "bond_angstrom");
        cJSON *ne = cJSON_GetObjectItem(j, "active_electrons");
        cJSON *no = cJSON_GetObjectItem(j, "active_orbitals");
        int ok = cJSON_IsArray(th) && cJSON_GetArraySize(th) == u->nexc && cJSON_IsNumber(b) &&
                 fabs(b->valuedouble - cfg->bond) < 1e-12 && cJSON_IsNumber(ne) && ne->valueint == cfg->nel &&
                 cJSON_IsNumber(no) && no->valueint == cfg->norb;
        if (ok) {
            for (int k = 0; k < u->nexc; k++) theta[k] = cJSON_GetArrayItem(th, k)->valuedouble;
            StateVec sv;
            sv_init(&sv, u->nq);
            *e_vqe = vqe_energy(u, ham, theta, &sv);
            sv_free(&sv);
            printf("  VQE parameters   : cached in %s\n", path);
            cJSON_Delete(j);
            return theta;
        }
        cJSON_Delete(j);
    }
    printf("  VQE parameters   : optimising UCCSD with COBYLA on the simulator (about 30 s)...\n");
    fflush(stdout);
    VQEOptions vo = {4000, 0.1, 1e-9, 0};
    VQEResult vr;
    vqe_optimize_local(u, ham, &vo, &vr);
    memcpy(theta, vr.theta, sizeof(double) * (size_t)u->nexc);
    *e_vqe = vr.energy;
    cJSON *j = cJSON_CreateObject();
    cJSON_AddNumberToObject(j, "bond_angstrom", cfg->bond);
    cJSON_AddNumberToObject(j, "active_electrons", cfg->nel);
    cJSON_AddNumberToObject(j, "active_orbitals", cfg->norb);
    cJSON_AddNumberToObject(j, "vqe_energy_Eh", vr.energy);
    cJSON_AddItemToObject(j, "theta", cJSON_CreateDoubleArray(vr.theta, vr.nparam));
    char *s = cJSON_Print(j);
    qc_mkdir_p(EXP_DIR);
    results_write_json(path, s);
    free(s);
    cJSON_Delete(j);
    vqe_result_free(&vr);
    return theta;
}

typedef struct {
    ExpConfig cfg;
    CasSystem cs;
    double e_rhf, e_fci, e_vqe, vqe_fidelity;
    UCCSD u;
    double *theta;
    ExpArm arms[3];
    int narms;
    IBMCalibration cal;
    IBMLayout lay;
    Circuit native[3];
    IBMLayout arm_lay[3];
    int offline;
} Plan;

static int parse_common(int argc, char **argv, char *cfgpath, size_t n, const char **calib, int *yes)
{
    snprintf(cfgpath, n, "%s", EXP_DEFAULT_CONFIG);
    calib[0] = calib[1] = NULL;
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--calibration") == 0 && i + 2 < argc) { calib[0] = argv[++i]; calib[1] = argv[++i]; }
        else if (strcmp(argv[i], "--yes") == 0) { if (yes) *yes = 1; }
        else if (argv[i][0] != '-') snprintf(cfgpath, n, "%s", argv[i]);
        else { fprintf(stderr, "error: unknown option '%s'\n", argv[i]); return -1; }
    }
    return 0;
}

static int build_plan(Plan *P, const char *cfgpath, const char *const *calib, IBMSession *sess)
{
    char err[512];
    if (exp_config_load(cfgpath, &P->cfg, err, sizeof err) != 0) { fprintf(stderr, "error: %s\n", err); return -1; }
    const ExpConfig *c = &P->cfg;
    printf("== experiment %s  (config %s)\n", c->name, cfgpath);
    cas_system_load(&P->cs, c->nel, c->norb, c->bond);
    FCIResult f;
    fci_solve(&f, &P->cs.fh, &P->cs.as, c->nel);
    P->e_fci = f.ground;
    fci_free(&f);
    P->e_rhf = P->cs.rhf.energy;
    uccsd_build(&P->u, P->cs.fh.nspin, c->nel);
    P->theta = exp_theta(c, &P->u, &P->cs.fh.ham, &P->e_vqe);
    exp_build_arms(c, &P->u, P->theta, &P->cs.fh.ham, P->arms, &P->narms, &P->vqe_fidelity);

    if (calib[0]) {
        size_t n1, n2;
        char *cf = qc_read_file(calib[0], &n1), *pr = qc_read_file(calib[1], &n2);
        if (!cf || !pr || ibm_calibration_parse(cf, pr, &P->cal) != 0) {
            fprintf(stderr, "error: cannot read calibration files\n");
            return -1;
        }
        free(cf);
        free(pr);
        if (layout_select(&P->cal, P->cs.fh.nspin, &P->lay) != 0) { fprintf(stderr, "error: no layout\n"); return -1; }
        P->offline = 1;
    } else {
        if (hw_get_session(sess) != 0) return -1;
        if (hw_select_backend(sess, c->backend, &P->cal, &P->lay, 1) != 0) return -1;
    }
    /* every arm on the same physical qubits with the same initial layout */
    for (int i = 0; i < P->narms; i++) {
        P->arm_lay[i] = P->lay;
        layout_compile(&P->arms[i].logical, &P->cal, &P->arm_lay[i], &P->native[i]);
    }
    return 0;
}

static void print_plan(const Plan *P)
{
    const ExpConfig *c = &P->cfg;
    printf("  system           : H2, R = %.3f A, cc-pVQZ; CAS(%d,%d) -> %d qubits, %zu Pauli terms\n",
           c->bond, c->nel, c->norb, P->cs.fh.nspin, P->cs.fh.ham.n);
    printf("  references       : RHF (full basis) %.8f | CAS-FCI %.8f | local VQE %.8f Eh\n", P->e_rhf, P->e_fci, P->e_vqe);
    printf("  backend          : %s%s, physical qubits [", P->cal.backend_name, P->offline ? " (offline calibration file)" : "");
    for (int l = 0; l < P->lay.nl; l++) printf("%d%s", P->lay.phys[l], l + 1 < P->lay.nl ? " " : "");
    printf("]\n");
    printf("  estimator        : one job, %d PUBs; precision %.3g, resilience_level %d, seed %d\n",
           P->narms, c->precision, c->resilience, c->seed);
    if (!c->ansatz_uccsd)
        printf("  treatment circuit: Givens, truncation %.0e, fidelity vs VQE state %.8f\n", c->gtrunc, P->vqe_fidelity);
    printf("\n  %-11s %-36s %5s %6s %8s %8s %15s\n", "arm", "purpose", "CZ", "depth", "runtime", "success", "exact E (Eh)");
    for (int i = 0; i < P->narms; i++) {
        double fid, dur, t2;
        hw_quality_estimate(&P->cal, &P->native[i], &P->arm_lay[i], &fid, &dur, &t2);
        printf("  %-11s %-36s %5d %6d %6.1fus %8.2f %15.8f\n", P->arms[i].name, P->arms[i].purpose,
               circ_count_kind(&P->native[i], G_CZ), circ_depth(&P->native[i]), dur, fid, P->arms[i].exact);
    }
    printf("\n  success = product of (1 - error) over CZ gates and readout: a heuristic, not a measurement\n");
}

static cJSON *plan_json(const Plan *P)
{
    const ExpConfig *c = &P->cfg;
    cJSON *j = cJSON_CreateObject();
    cJSON_AddStringToObject(j, "experiment", c->name);
    cJSON_AddStringToObject(j, "config_file", c->path);
    cJSON *cf = cJSON_AddObjectToObject(j, "config");
    cJSON_AddNumberToObject(cf, "bond_angstrom", c->bond);
    cJSON_AddNumberToObject(cf, "active_electrons", c->nel);
    cJSON_AddNumberToObject(cf, "active_orbitals", c->norb);
    cJSON_AddStringToObject(cf, "ansatz", c->ansatz_uccsd ? "uccsd" : "givens");
    cJSON_AddNumberToObject(cf, "givens_trunc", c->gtrunc);
    cJSON_AddNumberToObject(cf, "precision", c->precision);
    cJSON_AddNumberToObject(cf, "resilience_level", c->resilience);
    cJSON_AddNumberToObject(cf, "seed", c->seed);
    cJSON_AddStringToObject(j, "backend", P->cal.backend_name);
    cJSON_AddStringToObject(j, "backend_version", P->cal.backend_version);
    cJSON_AddItemToObject(j, "physical_qubits", cJSON_CreateIntArray(P->lay.phys, P->lay.nl));
    cJSON *ref = cJSON_AddObjectToObject(j, "reference");
    cJSON_AddNumberToObject(ref, "rhf_Eh", P->e_rhf);
    cJSON_AddNumberToObject(ref, "fci_Eh", P->e_fci);
    cJSON_AddNumberToObject(ref, "vqe_local_Eh", P->e_vqe);
    cJSON_AddNumberToObject(ref, "treatment_fidelity", P->vqe_fidelity);
    cJSON *arms = cJSON_AddArrayToObject(j, "arms");
    for (int i = 0; i < P->narms; i++) {
        cJSON *a = cJSON_CreateObject();
        cJSON_AddStringToObject(a, "name", P->arms[i].name);
        cJSON_AddStringToObject(a, "purpose", P->arms[i].purpose);
        cJSON_AddNumberToObject(a, "exact_Eh", P->arms[i].exact);
        cJSON_AddNumberToObject(a, "native_cz", circ_count_kind(&P->native[i], G_CZ));
        cJSON_AddNumberToObject(a, "native_depth", circ_depth(&P->native[i]));
        cJSON_AddItemToObject(a, "final_physical_positions", cJSON_CreateIntArray(P->arm_lay[i].final_pos, P->arm_lay[i].nl));
        cJSON_AddItemToArray(arms, a);
    }
    char *meta = results_metadata_json();
    cJSON_AddItemToObject(j, "metadata", cJSON_Parse(meta));
    free(meta);
    return j;
}

static int cmd_plan(int argc, char **argv, int run)
{
    char cfgpath[512];
    const char *calib[2];
    int yes = 0;
    if (parse_common(argc, argv, cfgpath, sizeof cfgpath, calib, &yes) != 0) return 1;
    Plan *P = qc_calloc(1, sizeof *P);
    IBMSession s;
    memset(&s, 0, sizeof s);
    if (build_plan(P, cfgpath, calib, &s) != 0) return 1;
    print_plan(P);

    qc_mkdir_p(EXP_DIR);
    cJSON *pj = plan_json(P);
    char path[600];
    snprintf(path, sizeof path, EXP_DIR "/%s_plan.json", P->cfg.name);
    char *txt = cJSON_Print(pj);
    results_write_json(path, txt);
    free(txt);
    if (!run) {
        printf("\nPlan only: nothing was submitted. Submit with: qchem-vqe experiment run %s\n", cfgpath);
        cJSON_Delete(pj);
        return 0;
    }
    if (P->offline) { fprintf(stderr, "error: an offline calibration file cannot be used to submit\n"); return 1; }

    int nobs = (int)P->cs.fh.ham.n;
    int ncz = 0;
    for (int i = 0; i < P->narms; i++) ncz += circ_count_kind(&P->native[i], G_CZ);
    if (!yes && !hw_confirm_submit(P->cal.backend_name, nobs * P->narms, ncz)) {
        printf("aborted: no job submitted\n");
        cJSON_Delete(pj);
        return 0;
    }
    const char *qasms[3];
    char **labels[3];
    for (int i = 0; i < P->narms; i++) {
        qasms[i] = qasm3_string(&P->native[i], P->cal.nqubit, 1);
        labels[i] = hw_pauli_labels(&P->cs.fh.ham, &P->arm_lay[i], P->cal.nqubit);
    }
    char *resp = NULL;
    int rc = ibm_job_submit_pubs(&s, P->cal.backend_name, qasms, P->narms, (const char *const *const *)labels, nobs,
                                 P->cfg.precision, P->cfg.resilience, P->cfg.seed, &resp);
    if (rc != 0) return 1;
    cJSON *rj = cJSON_Parse(resp);
    cJSON *jid = cJSON_GetObjectItem(rj, "id");
    if (!cJSON_IsString(jid)) { fprintf(stderr, "error: no job id in response\n"); return 1; }
    char job_id[128];
    snprintf(job_id, sizeof job_id, "%s", jid->valuestring);
    cJSON_Delete(rj);
    free(resp);

    /* record: everything needed to analyse the job later */
    cJSON_AddStringToObject(pj, "job_id", job_id);
    cJSON *obs = cJSON_AddArrayToObject(pj, "observable_coefficients");
    for (int t = 0; t < nobs; t++) cJSON_AddItemToArray(obs, cJSON_CreateNumber(creal(P->cs.fh.ham.t[t].c)));
    cJSON *arms = cJSON_GetObjectItem(pj, "arms");
    int i = 0;
    for (cJSON *a = arms->child; a; a = a->next, i++) {
        cJSON *lab = cJSON_AddArrayToObject(a, "observables");
        for (int t = 0; t < nobs; t++) cJSON_AddItemToArray(lab, cJSON_CreateString(labels[i][t]));
        char qp[600];
        snprintf(qp, sizeof qp, EXP_DIR "/%s_%s.qasm", job_id, P->arms[i].name);
        results_write_json(qp, qasms[i]);
    }
    snprintf(path, sizeof path, EXP_DIR "/%s.json", job_id);
    txt = cJSON_Print(pj);
    results_write_json(path, txt);
    free(txt);
    cJSON_Delete(pj);
    printf("\nSubmitted job %s (%d PUBs) to %s\n", job_id, P->narms, P->cal.backend_name);
    printf("Results, once completed:  qchem-vqe experiment results %s\n", job_id);
    return 0;
}

/* --- results --- */

static int arm_energy(cJSON *item, cJSON *coeffs, double *E, double *sigma)
{
    cJSON *res = cJSON_GetObjectItem(item, "results");
    cJSON *evs = res ? cJSON_GetObjectItem(res, "evs") : NULL;
    cJSON *stds = res ? cJSON_GetObjectItem(res, "stds") : NULL;
    cJSON *ed = evs ? cJSON_GetObjectItem(evs, "data") : NULL;
    cJSON *sd = stds ? cJSON_GetObjectItem(stds, "data") : NULL;
    if (!cJSON_IsString(ed) || !cJSON_IsString(sd)) return -1;
    double *ev = NULL, *st = NULL;
    size_t ne = 0, ns = 0;
    if (base64_decode_f64(ed->valuestring, &ev, &ne) != 0 || base64_decode_f64(sd->valuestring, &st, &ns) != 0) return -1;
    int nobs = cJSON_GetArraySize(coeffs);
    if ((int)ne != nobs || (int)ns != nobs) { free(ev); free(st); return -1; }
    double e = 0.0, v = 0.0;
    int k = 0;
    for (cJSON *c = coeffs->child; c; c = c->next, k++) {
        e += c->valuedouble * ev[k];
        v += c->valuedouble * c->valuedouble * st[k] * st[k];
    }
    free(ev);
    free(st);
    *E = e;
    *sigma = sqrt(v);
    return 0;
}

static int cmd_results(int argc, char **argv)
{
    if (argc < 1) { fprintf(stderr, "usage: qchem-vqe experiment results JOB_ID\n"); return 1; }
    const char *job = argv[0];
    char path[600];
    snprintf(path, sizeof path, EXP_DIR "/%s.json", job);
    size_t len = 0;
    char *txt = qc_read_file(path, &len);
    if (!txt) { fprintf(stderr, "error: no experiment record %s\n", path); return 1; }
    cJSON *rec = cJSON_Parse(txt);
    free(txt);

    IBMSession s;
    if (hw_get_session(&s) != 0) return 1;
    char *js = NULL;
    if (ibm_job_get(&s, job, &js) != 0) return 1;
    cJSON *jj = cJSON_Parse(js);
    free(js);
    cJSON *state = cJSON_GetObjectItem(jj, "state");
    cJSON *st = state ? cJSON_GetObjectItem(state, "status") : NULL;
    const char *status = cJSON_IsString(st) ? st->valuestring : "unknown";
    printf("== experiment %s, job %s: %s\n", cJSON_GetObjectItem(rec, "experiment")->valuestring, job, status);
    if (strcmp(status, "Completed") != 0) {
        cJSON *rs = state ? cJSON_GetObjectItem(state, "reason") : NULL;
        if (cJSON_IsString(rs)) printf("  reason: %s\n", rs->valuestring);
        cJSON_Delete(jj);
        cJSON_Delete(rec);
        return 0;
    }
    cJSON_Delete(jj);
    char *res = NULL;
    if (ibm_job_results(&s, job, &res) != 0) return 1;
    cJSON *rj = cJSON_Parse(res);
    free(res);
    cJSON *data = cJSON_GetObjectItem(rj, "data");
    cJSON *coeffs = cJSON_GetObjectItem(rec, "observable_coefficients");
    cJSON *arms = cJSON_GetObjectItem(rec, "arms");
    cJSON *ref = cJSON_GetObjectItem(rec, "reference");
    double e_fci = cJSON_GetObjectItem(ref, "fci_Eh")->valuedouble;

    printf("  backend %s, physical qubits as planned; one job, shared calibration window\n\n",
           cJSON_GetObjectItem(rec, "backend")->valuestring);
    printf("  %-11s %15s %24s %14s\n", "arm", "exact (Eh)", "QPU (Eh)", "QPU - exact");
    double E[3] = {0}, S[3] = {0}, X[3] = {0};
    const char *names[3] = {0};
    int n = 0;
    cJSON *out = cJSON_CreateObject();
    cJSON_AddStringToObject(out, "job_id", job);
    cJSON *oa = cJSON_AddArrayToObject(out, "arms");
    for (cJSON *a = arms->child; a && n < 3; a = a->next, n++) {
        names[n] = cJSON_GetObjectItem(a, "name")->valuestring;
        X[n] = cJSON_GetObjectItem(a, "exact_Eh")->valuedouble;
        if (arm_energy(cJSON_GetArrayItem(data, n), coeffs, &E[n], &S[n]) != 0) {
            fprintf(stderr, "error: unexpected result layout for PUB %d\n", n);
            return 1;
        }
        printf("  %-11s %15.8f %14.8f +/- %.6f %+14.6f\n", names[n], X[n], E[n], S[n], E[n] - X[n]);
        cJSON *o = cJSON_CreateObject();
        cJSON_AddStringToObject(o, "name", names[n]);
        cJSON_AddNumberToObject(o, "exact_Eh", X[n]);
        cJSON_AddNumberToObject(o, "qpu_Eh", E[n]);
        cJSON_AddNumberToObject(o, "qpu_sigma_Eh", S[n]);
        cJSON_AddItemToArray(oa, o);
    }
    int iv = -1, im = -1;
    for (int i = 0; i < n; i++) {
        if (strcmp(names[i], "vqe") == 0) iv = i;
        if (strcmp(names[i], "matched_hf") == 0) im = i;
    }
    if (iv >= 0 && im >= 0) {
        double dq = E[iv] - E[im], dx = X[iv] - X[im];
        double sd = sqrt(S[iv] * S[iv] + S[im] * S[im]);
        double corrected = X[im] + dq;
        printf("\n  correlation energy  E(vqe) - E(matched_hf):  exact %+.6f  QPU %+.6f +/- %.6f Eh\n", dx, dq, sd);
        printf("  control-corrected VQE energy  E_HF + [E(vqe) - E(matched_hf)]_QPU = %.8f +/- %.6f Eh\n", corrected, sd);
        printf("                                (CAS-FCI %.8f Eh, difference %+.6f Eh)\n", e_fci, corrected - e_fci);
        printf("  assumes the noise bias of the treatment equals that of its angle-zero control\n");
        cJSON_AddNumberToObject(out, "correlation_exact_Eh", dx);
        cJSON_AddNumberToObject(out, "correlation_qpu_Eh", dq);
        cJSON_AddNumberToObject(out, "correlation_qpu_sigma_Eh", sd);
        cJSON_AddNumberToObject(out, "control_corrected_vqe_Eh", corrected);
    }
    cJSON_AddItemToObject(out, "record", cJSON_Duplicate(rec, 1));
    snprintf(path, sizeof path, EXP_DIR "/%s_results.json", job);
    char *t2 = cJSON_Print(out);
    results_write_json(path, t2);
    free(t2);
    cJSON_Delete(out);
    cJSON_Delete(rj);
    cJSON_Delete(rec);
    return 0;
}

int cmd_experiment(int argc, char **argv)
{
    if (argc < 1) {
        fprintf(stderr, "usage:\n"
                        "  qchem-vqe experiment plan    [CONFIG] [--calibration CONF.json PROPS.json]\n"
                        "  qchem-vqe experiment run     [CONFIG] [--yes]\n"
                        "  qchem-vqe experiment results JOB_ID\n"
                        "CONFIG defaults to " EXP_DEFAULT_CONFIG "\n");
        return 1;
    }
    if (strcmp(argv[0], "plan") == 0) return cmd_plan(argc - 1, argv + 1, 0);
    if (strcmp(argv[0], "run") == 0) return cmd_plan(argc - 1, argv + 1, 1);
    if (strcmp(argv[0], "results") == 0) return cmd_results(argc - 1, argv + 1);
    fprintf(stderr, "unknown experiment subcommand '%s'\n", argv[0]);
    return 1;
}
