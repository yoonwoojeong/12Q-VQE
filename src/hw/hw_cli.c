/*
 * Hardware commands: hardware list|inspect|preflight, job status|results,
 * and vqe --mode hardware [--submit].
 *
 * A real QPU job is sent only with --submit AND an interactive confirmation
 * (type SUBMIT), or with --submit --yes.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cjson/cJSON.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define HW_DIR "results/hardware"
#define NLOG 12

typedef struct {
    char backend[128];
    char results_dir[512];
    double precision;
    double bond;
    int candidates;
    char job_id[128];
    int yes;
    int ansatz_uccsd;   /* 1: deep UCCSD circuit; 0 (default): shallow Givens preparation */
    double gtrunc;      /* natural-orbital weight that may be dropped by the Givens circuit */
} HwOpts;

static void hw_opts_parse(int argc, char **argv, HwOpts *o)
{
    memset(o, 0, sizeof *o);
    qc_results_dir(o->results_dir, sizeof o->results_dir);
    o->precision = 0.02;
    o->bond = 0.735;
    o->candidates = 8;
    o->gtrunc = 1e-4;
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--ansatz") == 0 && i + 1 < argc) {
            const char *a = argv[++i];
            if (strcmp(a, "uccsd") == 0) o->ansatz_uccsd = 1;
            else if (strcmp(a, "givens") == 0) o->ansatz_uccsd = 0;
            else qc_die("--ansatz must be 'givens' or 'uccsd'");
            continue;
        }
        if (strcmp(argv[i], "--givens-trunc") == 0 && i + 1 < argc) { o->gtrunc = atof(argv[++i]); continue; }
        if (strcmp(argv[i], "--results") == 0 && i + 1 < argc) { snprintf(o->results_dir, sizeof o->results_dir, "%s", argv[++i]); continue; }
        if (strcmp(argv[i], "--backend") == 0 && i + 1 < argc) snprintf(o->backend, sizeof o->backend, "%s", argv[++i]);
        else if (strcmp(argv[i], "--precision") == 0 && i + 1 < argc) o->precision = atof(argv[++i]);
        else if (strcmp(argv[i], "--yes") == 0) o->yes = 1;
        else if (argv[i][0] != '-' && o->job_id[0] == 0) snprintf(o->job_id, sizeof o->job_id, "%s", argv[i]);
    }
}

int hw_get_session(IBMSession *s)
{
    if (ibm_session_env(s) != 0) return -1;
    if (ibm_session_authenticate(s) != 0) {
        fprintf(stderr, "error: IBM IAM authentication failed (check the key with `qchem-vqe auth status`)\n");
        return -1;
    }
    return 0;
}

static void ensure_dir(const HwOpts *o)
{
    char d[1024];
    snprintf(d, sizeof d, "%s/hardware", o->results_dir);
    qc_mkdir_p(d);
}

/* Select the backend: explicit --backend, or the best online QPU by layout score. */
static int select_backend(IBMSession *s, const HwOpts *o, IBMCalibration *cal, IBMLayout *lay, int verbose)
{
    if (o->backend[0]) {
        char *cfg = NULL, *props = NULL, *st = NULL;
        if (ibm_get_backend_raw(s, o->backend, &cfg, &props, &st) != 0) qc_die("cannot read backend '%s'", o->backend);
        if (ibm_calibration_parse(cfg, props, cal) != 0) qc_die("cannot parse calibration of '%s'", o->backend);
        if (layout_select(cal, NLOG, lay) != 0) qc_die("no connected %d-qubit layout on '%s'", NLOG, o->backend);
        free(cfg); free(props); free(st);
        return 0;
    }
    IBMBackendInfo *bi = NULL;
    int nb = 0;
    if (ibm_list_backends(s, &bi, &nb) != 0) return -1;
    int best_idx = -1;
    double best = 1e300;
    IBMCalibration best_cal;
    IBMLayout best_lay;
    int tried = 0;
    for (int i = 0; i < nb && tried < o->candidates; i++) {
        if (bi[i].is_simulator || strcmp(bi[i].status, "online") != 0 || bi[i].n_qubits < NLOG) continue;
        tried++;
        char *cfg = NULL, *props = NULL, *st = NULL;
        if (ibm_get_backend_raw(s, bi[i].name, &cfg, &props, &st) != 0) { free(cfg); free(props); free(st); continue; }
        IBMCalibration c;
        IBMLayout l;
        int ok = ibm_calibration_parse(cfg, props, &c) == 0 && layout_select(&c, NLOG, &l) == 0;
        if (verbose && ok) printf("  candidate %-18s %3d qubits  layout error score %.5f  queue %d\n",
                                  bi[i].name, bi[i].n_qubits, l.score, bi[i].queue_length);
        if (ok && l.score < best) {
            best = l.score;
            best_idx = i;
            best_cal = c;
            best_lay = l;
        }
        free(cfg); free(props); free(st);
    }
    if (best_idx < 0) {
        free(bi);
        fprintf(stderr, "error: no online QPU with a connected %d-qubit layout and CZ gates\n", NLOG);
        return -1;
    }
    *cal = best_cal;
    *lay = best_lay;
    printf("  selected backend  : %s (%d qubits, queue %d, layout error score %.5f)\n",
           bi[best_idx].name, bi[best_idx].n_qubits, bi[best_idx].queue_length, best);
    free(bi);
    return 0;
}

static void print_layout(const IBMCalibration *cal, const IBMLayout *lay)
{
    printf("  physical layout   : logical 0..11 -> physical [");
    for (int l = 0; l < NLOG; l++) printf("%d%s", lay->phys[l], l + 1 < NLOG ? ", " : "");
    printf("]\n");
    printf("  layout error score: %.5f (sum readout + CZ errors on the subgraph)\n", lay->score);
    (void)cal;
}

/* Heuristic quality of a compiled circuit: product of (1 - error) over CZ gates and
   readout of the layout qubits, and runtime as depth x mean CZ length. Not a measurement. */
void hw_quality_estimate(const IBMCalibration *cal, const Circuit *native, const IBMLayout *lay,
                             double *fid, double *dur_us, double *min_t2_us)
{
    double f = 1.0, len_s = 0.0, len_n = 0.0;
    for (int e = 0; e < cal->nedge; e++)
        if (!isnan(cal->edge_len_ns[e])) { len_s += cal->edge_len_ns[e]; len_n++; }
    for (int k = 0; k < native->ng; k++) {
        const Gate *g = &native->g[k];
        if (g->kind != G_CZ) continue;
        for (int e = 0; e < cal->nedge; e++)
            if ((cal->edge_a[e] == g->q0 && cal->edge_b[e] == g->q1) || (cal->edge_a[e] == g->q1 && cal->edge_b[e] == g->q0)) {
                f *= 1.0 - (isnan(cal->edge_err[e]) ? 0.01 : cal->edge_err[e]);
                break;
            }
    }
    double t2 = 1e300;
    for (int i = 0; i < lay->nsubset; i++) {
        int q = lay->subset[i];
        f *= 1.0 - (isnan(cal->readout_err[q]) ? 0.02 : cal->readout_err[q]);
        if (!isnan(cal->t2_us[q]) && cal->t2_us[q] < t2) t2 = cal->t2_us[q];
    }
    *fid = f;
    (void)len_n;
    (void)len_s;
    *dur_us = circ_duration_ns(native, cal) / 1000.0;
    *min_t2_us = t2 < 1e300 ? t2 : NAN;
}

/* Logical circuit that prepares the optimised VQE state.
   givens (default): natural-orbital Givens preparation of the UCCSD state, verified on the simulator.
   uccsd: the UCCSD circuit itself. *e_circ is the noiseless energy of the returned circuit. */
static void build_logical_circuit(const HwOpts *o, const UCCSD *u, const double *theta, const PauliSum *ham,
                                  Circuit *out, double *e_circ, double *fidelity)
{
    StateVec target, v;
    sv_init(&target, u->nq);
    sv_init(&v, u->nq);
    uccsd_prepare(u, theta, &target);
    if (o->ansatz_uccsd) {
        uccsd_circuit(u, theta, out);
        *fidelity = 1.0;
        printf("  ansatz            : UCCSD circuit (%d excitations)\n", u->nexc);
    } else {
        GivensInfo gi;
        givens_state_circuit(&target, u->nq, o->gtrunc, out, &gi);
        *fidelity = gi.fidelity;
        printf("  ansatz            : Givens preparation of the optimised UCCSD state\n");
        printf("                      %d of %d natural-orbital pairs (weight kept %.8f), %d pair + %d orbital rotations\n",
               gi.r, u->nq / 2, gi.kept_weight, gi.n_pair_rot, gi.n_givens);
        printf("                      fidelity with the VQE state on the simulator: %.8f\n", gi.fidelity);
        if (gi.fidelity < 0.999) qc_die("Givens circuit fidelity %.6f is below 0.999; lower --givens-trunc", gi.fidelity);
    }
    circ_simulate(out, &v);
    *e_circ = sv_expect_sum(&v, ham);
    printf("  noiseless energy  : %.10f Eh for this circuit (VQE state %.10f Eh)\n", *e_circ, sv_expect_sum(&target, ham));
    sv_free(&target);
    sv_free(&v);
}

/* ---------- hardware list / inspect / preflight ---------- */

static int cmd_list(const HwOpts *o)
{
    (void)o;
    IBMSession s;
    if (hw_get_session(&s) != 0) return 1;
    IBMBackendInfo *bi = NULL;
    int nb = 0;
    if (ibm_list_backends(&s, &bi, &nb) != 0) return 1;
    printf("%-22s %-10s %-8s %-7s %s\n", "backend", "status", "qubits", "queue", "simulator");
    for (int i = 0; i < nb; i++)
        printf("%-22s %-10s %-8d %-7d %s\n", bi[i].name, bi[i].status, bi[i].n_qubits, bi[i].queue_length,
               bi[i].is_simulator ? "yes" : "no");
    free(bi);
    return 0;
}

static void save_calibration_snapshot(const HwOpts *o, const char *backend, const char *cfg, const char *props)
{
    ensure_dir(o);
    char p[1024];
    snprintf(p, sizeof p, "%s/hardware/%s_configuration.json", o->results_dir, backend);
    results_write_json(p, cfg);
    snprintf(p, sizeof p, "%s/hardware/%s_properties.json", o->results_dir, backend);
    results_write_json(p, props);
}

static int cmd_inspect(const HwOpts *o, int preflight)
{
    IBMSession s;
    if (hw_get_session(&s) != 0) return 1;
    IBMCalibration cal;
    IBMLayout lay;
    if (select_backend(&s, o, &cal, &lay, 1) != 0) return 1;

    /* re-fetch raw documents for the snapshot */
    char *cfg = NULL, *props = NULL, *st = NULL;
    ibm_get_backend_raw(&s, cal.backend_name, &cfg, &props, &st);
    save_calibration_snapshot(o, cal.backend_name, cfg, props);

    printf("== backend %s (version %s)\n", cal.backend_name, cal.backend_version);
    printf("  qubits %d, couplers %d, basis: %s%s\n", cal.nqubit, cal.nedge, cal.has_cz ? "cz " : "", cal.has_ecr ? "ecr" : "");
    if (st) {
        cJSON *j = cJSON_Parse(st);
        cJSON *sv = cJSON_GetObjectItem(j, "status");
        cJSON *q = cJSON_GetObjectItem(j, "length_queue");
        printf("  status: %s, queue length %d\n", cJSON_IsString(sv) ? sv->valuestring : "?", cJSON_IsNumber(q) ? q->valueint : -1);
        cJSON_Delete(j);
    }
    /* calibration statistics */
    double cz_s = 0, cz_n = 0, cz_max = 0, ro_s = 0, ro_n = 0, len_s = 0, len_n = 0, sx_s = 0, sx_n = 0;
    for (int e = 0; e < cal.nedge; e++) {
        if (!isnan(cal.edge_err[e])) { cz_s += cal.edge_err[e]; cz_n++; if (cal.edge_err[e] > cz_max) cz_max = cal.edge_err[e]; }
        if (!isnan(cal.edge_len_ns[e])) { len_s += cal.edge_len_ns[e]; len_n++; }
    }
    for (int q = 0; q < cal.nqubit; q++) {
        if (!isnan(cal.readout_err[q])) { ro_s += cal.readout_err[q]; ro_n++; }
        if (!isnan(cal.sx_err[q])) { sx_s += cal.sx_err[q]; sx_n++; }
    }
    printf("  two-qubit (CZ) error: mean %.4f  max %.4f  (%0.f couplers)\n", cz_n ? cz_s / cz_n : NAN, cz_max, cz_n);
    printf("  CZ gate length      : mean %.1f ns (%0.f couplers)\n", len_n ? len_s / len_n : NAN, len_n);
    printf("  readout error       : mean %.4f (%0.f qubits)\n", ro_n ? ro_s / ro_n : NAN, ro_n);
    printf("  SX error            : mean %.5f (%0.f qubits)\n", sx_n ? sx_s / sx_n : NAN, sx_n);
    print_layout(&cal, &lay);
    printf("  layout qubits' T1/T2 (us): ");
    for (int i = 0; i < lay.nsubset && i < 12; i++) {
        int q = lay.subset[i];
        printf("q%d %.0f/%.0f  ", q, cal.t1_us[q], cal.t2_us[q]);
    }
    printf("\n");
    if (preflight) {
        /* circuit cost on this layout, from the optimised local VQE parameters */
        CasSystem cs;
        cas_system_load(&cs, 2, 6, o->bond);
        UCCSD u;
        uccsd_build(&u, cs.fh.nspin, 2);
        FILE *f = NULL;
        char tpath[1024];
        snprintf(tpath, sizeof tpath, "%s/vqe_local_CAS2_6.json", o->results_dir);
        f = fopen(tpath, "rb");
        double *theta = qc_calloc((size_t)u.nexc, sizeof(double));
        if (f) {
            size_t len = 0;
            char *txt = qc_read_file(tpath, &len);
            cJSON *j = cJSON_Parse(txt);
            cJSON *th = cJSON_GetObjectItem(j, "theta");
            for (int k = 0; k < u.nexc && cJSON_IsArray(th) && k < cJSON_GetArraySize(th); k++)
                theta[k] = cJSON_GetArrayItem(th, k)->valuedouble;
            cJSON_Delete(j);
            free(txt);
            fclose(f);
        } else {
            printf("  (no results/vqe_local_CAS2_6.json: using theta=0 for the cost estimate)\n");
        }
        Circuit logical, native;
        double e_circ_pf, fid_pf;
        build_logical_circuit(o, &u, theta, &cs.fh.ham, &logical, &e_circ_pf, &fid_pf);
        IBMLayout l2 = lay;
        layout_compile(&logical, &cal, &l2, &native);
        int ncz = circ_count_kind(&native, G_CZ);
        double fid = 1.0;
        for (int k = 0; k < native.ng; k++) {
            const Gate *g = &native.g[k];
            if (g->kind != G_CZ) continue;
            for (int e = 0; e < cal.nedge; e++) {
                if ((cal.edge_a[e] == g->q0 && cal.edge_b[e] == g->q1) || (cal.edge_a[e] == g->q1 && cal.edge_b[e] == g->q0)) {
                    double err = isnan(cal.edge_err[e]) ? 0.01 : cal.edge_err[e];
                    fid *= 1.0 - err;
                    break;
                }
            }
        }
        for (int i = 0; i < l2.nsubset; i++) {
            double r = isnan(cal.readout_err[l2.subset[i]]) ? 0.02 : cal.readout_err[l2.subset[i]];
            fid *= 1.0 - r;
        }
        double tns = circ_duration_ns(&native, &cal);
        printf("== preflight (circuit on this layout, theta from results/vqe_local_CAS2_6.json)\n");
        printf("  logical gates      : %d (CX %d, depth %d)\n", logical.ng, circ_count_kind(&logical, G_CX), circ_depth(&logical));
        printf("  routing SWAPs      : %d (3 CZ each)\n", l2.nswaps);
        printf("  native gates       : %d, CZ %d, native depth %d\n", native.ng, ncz, circ_depth(&native));
        printf("  est. duration      : %.1f us (timed critical path: rz virtual, calibrated sx/cz lengths)\n", tns / 1000.0);
        double min_t1 = 1e9, min_t2 = 1e9;
        for (int i = 0; i < l2.nsubset; i++) {
            int q = l2.subset[i];
            if (!isnan(cal.t1_us[q]) && cal.t1_us[q] < min_t1) min_t1 = cal.t1_us[q];
            if (!isnan(cal.t2_us[q]) && cal.t2_us[q] < min_t2) min_t2 = cal.t2_us[q];
        }
        printf("  worst T1 / T2 in layout: %.0f / %.0f us\n", min_t1, min_t2);
        printf("  heuristic fidelity (product of 1-error): %.3e  -- circuit-level estimate, not a measurement\n", fid);
        printf("  Pauli observables  : %zu (estimator groups commuting terms internally)\n", cs.fh.ham.n);
        if (tns / 1000.0 > (min_t2 > 0 && !isnan(min_t2) ? min_t2 : 100.0))
            printf("  WARNING: estimated duration exceeds the worst T2 of the layout; expect decoherence-dominated results\n");
        free(theta);
        uccsd_free(&u);
        circ_free(&logical);
        circ_free(&native);
        cas_system_free(&cs);
    }
    free(cfg); free(props); free(st);
    return 0;
}

int cmd_hardware(int argc, char **argv)
{
    if (argc < 1) qc_die("usage: qchem-vqe hardware list | inspect | preflight [--backend NAME]");
    HwOpts o;
    hw_opts_parse(argc - 1, argv + 1, &o);
    if (strcmp(argv[0], "list") == 0) return cmd_list(&o);
    if (strcmp(argv[0], "inspect") == 0) return cmd_inspect(&o, 0);
    if (strcmp(argv[0], "preflight") == 0) return cmd_inspect(&o, 1);
    qc_die("unknown hardware subcommand '%s'", argv[0]);
}

/* ---------- job status / results ---------- */

static int cmd_job_status(const HwOpts *o)
{
    IBMSession s;
    if (hw_get_session(&s) != 0) return 1;
    char *js = NULL;
    if (ibm_job_get(&s, o->job_id, &js) != 0) qc_die("cannot read job %s", o->job_id);
    cJSON *j = cJSON_Parse(js);
    cJSON *state = cJSON_GetObjectItem(j, "state");
    cJSON *st = state ? cJSON_GetObjectItem(state, "status") : NULL;
    cJSON *rs = state ? cJSON_GetObjectItem(state, "reason") : NULL;
    cJSON *be = cJSON_GetObjectItem(j, "backend");
    printf("job %s\n  backend: %s\n  status : %s\n", o->job_id,
           cJSON_IsString(be) ? be->valuestring : "?", cJSON_IsString(st) ? st->valuestring : "?");
    if (cJSON_IsString(rs)) printf("  reason : %s\n", rs->valuestring);
    cJSON_Delete(j);
    free(js);
    return 0;
}

int ibm_parse_estimator_results_json(const char *res_json, double **ev_out, double **sd_out, size_t *nout)
{
    *ev_out = NULL;
    *sd_out = NULL;
    *nout = 0;
    cJSON *rj = cJSON_Parse(res_json);
    if (!rj) return -1;
    cJSON *data = cJSON_GetObjectItem(rj, "data");
    cJSON *item = cJSON_IsArray(data) ? cJSON_GetArrayItem(data, 0) : NULL;
    if (!item && cJSON_IsArray(rj)) item = cJSON_GetArrayItem(rj, 0);
    cJSON *results = item ? cJSON_GetObjectItem(item, "results") : NULL;
    if (!results && item) results = cJSON_GetObjectItem(item, "data");
    if (!results) results = cJSON_GetObjectItem(rj, "results");
    if (cJSON_IsArray(results)) results = cJSON_GetArrayItem(results, 0);
    if (results && cJSON_GetObjectItem(results, "data")) results = cJSON_GetObjectItem(results, "data");
    cJSON *evs = results ? cJSON_GetObjectItem(results, "evs") : NULL;
    cJSON *stds = results ? cJSON_GetObjectItem(results, "stds") : NULL;
    if (!evs || !stds) {
        cJSON_Delete(rj);
        return -2;
    }

    double *ev = NULL, *sd = NULL;
    size_t nev = 0, nsd = 0;
    if (cJSON_IsArray(evs) && cJSON_IsArray(stds)) {
        nev = (size_t)cJSON_GetArraySize(evs);
        nsd = (size_t)cJSON_GetArraySize(stds);
        ev = qc_calloc(nev ? nev : 1, sizeof(double));
        sd = qc_calloc(nsd ? nsd : 1, sizeof(double));
        int idx = 0;
        for (cJSON *e = evs->child; e; e = e->next, idx++) ev[idx] = e->valuedouble;
        idx = 0;
        for (cJSON *e = stds->child; e; e = e->next, idx++) sd[idx] = e->valuedouble;
    } else if (cJSON_IsObject(evs) && cJSON_GetObjectItem(evs, "data") &&
               cJSON_IsObject(stds) && cJSON_GetObjectItem(stds, "data")) {
        if (base64_decode_f64(cJSON_GetObjectItem(evs, "data")->valuestring, &ev, &nev) != 0 ||
            base64_decode_f64(cJSON_GetObjectItem(stds, "data")->valuestring, &sd, &nsd) != 0) {
            free(ev); free(sd); cJSON_Delete(rj); return -3;
        }
    } else {
        cJSON_Delete(rj);
        return -4;
    }
    cJSON_Delete(rj);
    if (nev != nsd) { free(ev); free(sd); return -5; }
    *ev_out = ev;
    *sd_out = sd;
    *nout = nev;
    return 0;
}

static int cmd_job_results(const HwOpts *o)
{
    IBMSession s;
    if (hw_get_session(&s) != 0) return 1;
    char *js = NULL;
    if (ibm_job_get(&s, o->job_id, &js) != 0) qc_die("cannot read job %s", o->job_id);
    cJSON *j = cJSON_Parse(js);
    cJSON *state = cJSON_GetObjectItem(j, "state");
    cJSON *st = state ? cJSON_GetObjectItem(state, "status") : NULL;
    const char *status = cJSON_IsString(st) ? st->valuestring : "unknown";
    printf("job %s status: %s\n", o->job_id, status);
    int completed = (strcmp(status, "Completed") == 0);
    cJSON_Delete(j);
    free(js);
    if (!completed) {
        printf("  results are available only for a Completed job\n");
        return 0;
    }

    /* submission record written at submit time */
    char rpath[1024];
    snprintf(rpath, sizeof rpath, "%s/hardware/%s.json", o->results_dir, o->job_id);
    size_t rlen = 0;
    char *rec_txt = qc_read_file(rpath, &rlen);
    if (!rec_txt) qc_die("no submission record %s (was this job submitted by this tool?)", rpath);
    cJSON *rec = cJSON_Parse(rec_txt);
    free(rec_txt);
    cJSON *obs = cJSON_GetObjectItem(rec, "observables");
    int nobs = cJSON_GetArraySize(obs);

    char *res = NULL;
    if (ibm_job_results(&s, o->job_id, &res) != 0) qc_die("cannot read results of %s", o->job_id);
    double *ev = NULL, *sd = NULL;
    size_t nev = 0;
    if (ibm_parse_estimator_results_json(res, &ev, &sd, &nev) != 0) {
        fprintf(stderr, "raw results JSON: %.600s\n", res);
        free(res);
        qc_die("unexpected result layout (no evs/stds)");
    }
    free(res);
    if ((int)nev != nobs) qc_die("result length %zu does not match %d observables", nev, nobs);

    double E = 0.0, var = 0.0;
    int k = 0;
    for (cJSON *o2 = obs->child; o2; o2 = o2->next, k++) {
        double c = cJSON_GetObjectItem(o2, "coeff")->valuedouble;
        E += c * ev[k];
        var += c * c * sd[k] * sd[k];
    }
    double sigma = sqrt(var);
    cJSON *ref = cJSON_GetObjectItem(rec, "reference");
    double fci = cJSON_GetObjectItem(ref, "fci_Eh")->valuedouble;
    double loc = cJSON_GetObjectItem(ref, "vqe_local_Eh")->valuedouble;
    double rhf = cJSON_GetObjectItem(ref, "rhf_Eh")->valuedouble;
    printf("  QPU energy (CAS(2,6)) : %.8f +/- %.8f Eh  (independent-term error propagation)\n", E, sigma);
    cJSON *nc = cJSON_GetObjectItem(ref, "noiseless_circuit_Eh");
    if (cJSON_IsNumber(nc))
        printf("  submitted circuit     : %.8f Eh noiseless  QPU - noiseless = %+.6f Eh\n", nc->valuedouble, E - nc->valuedouble);
    printf("  local VQE             : %.8f Eh   QPU - local = %+.6f Eh\n", loc, E - loc);
    printf("  CAS-FCI reference     : %.8f Eh   QPU - FCI   = %+.6f Eh\n", fci, E - fci);
    printf("  RHF (full basis)      : %.8f Eh\n", rhf);
    printf("  QPU correlation energy recovered: %.1f %% of FCI (vs RHF)\n",
           100.0 * (E - rhf) / (fci - rhf));

    cJSON *out = cJSON_CreateObject();
    cJSON_AddStringToObject(out, "job_id", o->job_id);
    cJSON_AddStringToObject(out, "scope", "CAS(2,6) active-space energy evaluated on IBM QPU");
    cJSON_AddNumberToObject(out, "qpu_energy_Eh", E);
    cJSON_AddNumberToObject(out, "qpu_energy_sigma_Eh", sigma);
    cJSON_AddNumberToObject(out, "vqe_local_Eh", loc);
    cJSON_AddNumberToObject(out, "fci_Eh", fci);
    cJSON_AddNumberToObject(out, "rhf_Eh", rhf);
    cJSON_AddNumberToObject(out, "n_observables", nobs);
    cJSON_AddItemToObject(out, "submission", cJSON_Duplicate(rec, 1));
    char *txt = cJSON_Print(out);
    snprintf(rpath, sizeof rpath, "%s/hardware/%s_results.json", o->results_dir, o->job_id);
    results_write_json(rpath, txt);
    free(txt);
    cJSON_Delete(out);
    cJSON_Delete(rec);
    free(ev);
    free(sd);
    return 0;
}

int cmd_job(int argc, char **argv)
{
    if (argc < 2) qc_die("usage: qchem-vqe job status JOB_ID | job results JOB_ID");
    HwOpts o;
    hw_opts_parse(argc - 1, argv + 1, &o);
    if (!o.job_id[0]) qc_die("missing JOB_ID");
    if (strcmp(argv[0], "status") == 0) return cmd_job_status(&o);
    if (strcmp(argv[0], "results") == 0) return cmd_job_results(&o);
    qc_die("unknown job subcommand '%s'", argv[0]);
}

/* ---------- vqe --mode hardware ---------- */

int hw_confirm_submit(const char *backend, int nobs, int ncz)
{
    printf("\nAbout to submit an Estimator V2 job to REAL hardware: backend %s,\n", backend);
    printf("%d observables, %d native CZ gates. This consumes IBM Quantum time.\n", nobs, ncz);
    printf("Type SUBMIT to continue: ");
    fflush(stdout);
    char line[64] = {0};
    if (!fgets(line, sizeof line, stdin)) return 0;
    line[strcspn(line, "\r\n")] = 0;
    return strcmp(line, "SUBMIT") == 0;
}

int hw_run_vqe_hardware(int argc, char **argv, int submit, int yes)
{
    HwOpts o;
    hw_opts_parse(argc, argv, &o);
    o.yes = yes;
    qc_print_banner("VQE on IBM Quantum hardware (evaluation of the optimised UCCSD state)");

    /* optimised local parameters */
    char tpath[1024];
    snprintf(tpath, sizeof tpath, "%s/vqe_local_CAS2_6.json", o.results_dir);
    size_t tl = 0;
    char *txt = qc_read_file(tpath, &tl);
    if (!txt) qc_die("missing %s: run `qchem-vqe vqe --mode local` first", tpath);
    cJSON *lj = cJSON_Parse(txt);
    free(txt);
    cJSON *th = cJSON_GetObjectItem(lj, "theta");
    double local_E = cJSON_GetObjectItem(lj, "vqe_local_Eh")->valuedouble;
    double fci_E = cJSON_GetObjectItem(lj, "fci_cas_Eh")->valuedouble;
    double rhf_E = cJSON_GetObjectItem(lj, "rhf_full_basis_Eh")->valuedouble;

    CasSystem cs;
    cas_system_load(&cs, 2, 6, o.bond);
    UCCSD u;
    uccsd_build(&u, cs.fh.nspin, 2);
    if (!cJSON_IsArray(th) || cJSON_GetArraySize(th) != u.nexc) qc_die("theta length does not match UCCSD (%d)", u.nexc);
    double *theta = qc_calloc((size_t)u.nexc, sizeof(double));
    for (int k = 0; k < u.nexc; k++) theta[k] = cJSON_GetArrayItem(th, k)->valuedouble;
    cJSON_Delete(lj);

    Circuit logical;
    double e_circ = 0.0, circ_fid = 1.0;
    build_logical_circuit(&o, &u, theta, &cs.fh.ham, &logical, &e_circ, &circ_fid);

    IBMSession s;
    if (hw_get_session(&s) != 0) return 1;
    IBMCalibration cal;
    IBMLayout lay;
    if (select_backend(&s, &o, &cal, &lay, 1) != 0) return 1;
    print_layout(&cal, &lay);

    Circuit native;
    layout_compile(&logical, &cal, &lay, &native);
    char *qasm = qasm3_string(&native, cal.nqubit, 1);

    /* observables: each JW Pauli term mapped to the physical register (qiskit label: qubit 0 rightmost) */
    int nobs = (int)cs.fh.ham.n;
    char **labels = qc_calloc((size_t)nobs, sizeof(char *));
    for (int t = 0; t < nobs; t++) {
        char *lab = qc_calloc((size_t)cal.nqubit + 1, 1);
        for (int p = 0; p < cal.nqubit; p++) lab[p] = 'I';
        for (int l = 0; l < NLOG; l++) {
            int xb = (cs.fh.ham.t[t].x >> l) & 1, zb = (cs.fh.ham.t[t].z >> l) & 1;
            char c = (xb && zb) ? 'Y' : (xb ? 'X' : (zb ? 'Z' : 0));
            if (!c) continue;
            int p = lay.final_pos[l];
            lab[cal.nqubit - 1 - p] = c;
        }
        lab[cal.nqubit] = 0;
        labels[t] = lab;
    }

    printf("  circuit           : %d logical gates -> %d native gates, CZ %d, native depth %d, SWAPs %d\n",
           logical.ng, native.ng, circ_count_kind(&native, G_CZ), circ_depth(&native), lay.nswaps);
    printf("  observables       : %d Pauli strings over %d physical qubits\n", nobs, cal.nqubit);
    printf("  estimator options : default_precision %.3f, resilience_level 1\n", o.precision);
    printf("  local VQE         : %.8f Eh  (FCI %.8f Eh)\n", local_E, fci_E);
    double q_fid, q_dur, q_t2;
    hw_quality_estimate(&cal, &native, &lay, &q_fid, &q_dur, &q_t2);
    printf("  est. success prob.: %.1e (product of 1-error over CZ and readout; heuristic)\n", q_fid);
    printf("  est. runtime      : %.0f us per shot vs worst T2 in layout %.0f us\n", q_dur, q_t2);
    if (q_fid < 1e-3 || (!isnan(q_t2) && q_dur > q_t2))
        printf("  WARNING: at this depth the QPU result is expected to be dominated by noise;\n"
               "           the measured energy will not reflect the VQE state.\n");

    if (!submit) {
        printf("\nDry run: nothing was submitted. Re-run with --submit to send the job after confirmation.\n");
        goto cleanup_ok;
    }
    if (!o.yes && !hw_confirm_submit(cal.backend_name, nobs, circ_count_kind(&native, G_CZ))) {
        printf("aborted: no job submitted\n");
        goto cleanup_ok;
    }
    char *resp = NULL;
    if (ibm_job_submit(&s, cal.backend_name, qasm, (const char *const *)labels, nobs, o.precision, &resp) != 0) {
        free(resp);
        goto cleanup_fail;
    }
    cJSON *rj = cJSON_Parse(resp);
    cJSON *jid = cJSON_GetObjectItem(rj, "id");
    if (!cJSON_IsString(jid)) { fprintf(stderr, "error: no job id in response\n"); cJSON_Delete(rj); free(resp); goto cleanup_fail; }
    char job_id[128];
    snprintf(job_id, sizeof job_id, "%s", jid->valuestring);
    cJSON_Delete(rj);
    free(resp);
    printf("  submitted job     : %s\n", job_id);

    /* submission record: everything needed to evaluate the job later */
    ensure_dir(&o);
    cJSON *rec = cJSON_CreateObject();
    cJSON_AddStringToObject(rec, "job_id", job_id);
    cJSON_AddStringToObject(rec, "backend", cal.backend_name);
    cJSON_AddStringToObject(rec, "backend_version", cal.backend_version);
    cJSON_AddNumberToObject(rec, "estimator_version", 2);
    cJSON_AddNumberToObject(rec, "resilience_level", 1);
    cJSON_AddNumberToObject(rec, "default_precision", o.precision);
    cJSON_AddNumberToObject(rec, "logical_gates", logical.ng);
    cJSON_AddNumberToObject(rec, "native_gates", native.ng);
    cJSON_AddNumberToObject(rec, "native_cz", circ_count_kind(&native, G_CZ));
    cJSON_AddNumberToObject(rec, "native_depth", circ_depth(&native));
    cJSON_AddNumberToObject(rec, "routing_swaps", lay.nswaps);
    cJSON *phys = cJSON_CreateIntArray(lay.phys, NLOG);
    cJSON_AddItemToObject(rec, "initial_physical_layout", phys);
    cJSON *fin = cJSON_CreateIntArray(lay.final_pos, NLOG);
    cJSON_AddItemToObject(rec, "final_physical_positions", fin);
    cJSON_AddNumberToObject(rec, "layout_error_score", lay.score);
    cJSON_AddStringToObject(rec, "ansatz", o.ansatz_uccsd ? "uccsd" : "givens");
    if (!o.ansatz_uccsd) cJSON_AddNumberToObject(rec, "givens_truncation", o.gtrunc);
    cJSON_AddNumberToObject(rec, "circuit_fidelity_vs_vqe_state", circ_fid);
    cJSON *ob = cJSON_AddArrayToObject(rec, "observables");
    for (int t = 0; t < nobs; t++) {
        cJSON *e = cJSON_CreateObject();
        cJSON_AddStringToObject(e, "label", labels[t]);
        cJSON_AddNumberToObject(e, "coeff", creal(cs.fh.ham.t[t].c));
        cJSON_AddItemToArray(ob, e);
    }
    cJSON *ref = cJSON_AddObjectToObject(rec, "reference");
    cJSON_AddNumberToObject(ref, "rhf_Eh", rhf_E);
    cJSON_AddNumberToObject(ref, "fci_Eh", fci_E);
    cJSON_AddNumberToObject(ref, "vqe_local_Eh", local_E);
    cJSON_AddNumberToObject(ref, "noiseless_circuit_Eh", e_circ);
    char *meta = results_metadata_json();
    cJSON_AddItemToObject(rec, "metadata", cJSON_Parse(meta));
    free(meta);
    char *rtxt = cJSON_Print(rec);
    char rpath[1024];
    snprintf(rpath, sizeof rpath, "%s/hardware/%s.json", o.results_dir, job_id);
    results_write_json(rpath, rtxt);
    snprintf(rpath, sizeof rpath, "%s/hardware/%s.qasm", o.results_dir, job_id);
    results_write_json(rpath, qasm);
    free(rtxt);
    cJSON_Delete(rec);
    printf("  check status with : qchem-vqe job status %s\n", job_id);
    printf("  fetch results with: qchem-vqe job results %s\n", job_id);
cleanup_ok:
    for (int t = 0; t < nobs; t++) free(labels[t]);
    free(labels);
    free(qasm);
    circ_free(&native);
    circ_free(&logical);
    free(theta);
    uccsd_free(&u);
    cas_system_free(&cs);
    return 0;
cleanup_fail:
    for (int t = 0; t < nobs; t++) free(labels[t]);
    free(labels);
    free(qasm);
    circ_free(&native);
    circ_free(&logical);
    free(theta);
    uccsd_free(&u);
    cas_system_free(&cs);
    return 1;
}

int hw_select_backend(IBMSession *s, const char *backend, IBMCalibration *cal, IBMLayout *lay, int verbose)
{
    HwOpts o;
    hw_opts_parse(0, NULL, &o);
    if (backend && *backend && strcmp(backend, "auto") != 0) snprintf(o.backend, sizeof o.backend, "%s", backend);
    return select_backend(s, &o, cal, lay, verbose);
}

/* Qiskit-ordered labels (qubit 0 rightmost) of every Hamiltonian term on the physical
   register, with logical qubit l at physical lay->final_pos[l]. Caller frees. */
char **hw_pauli_labels(const PauliSum *ham, const IBMLayout *lay, int nphys)
{
    char **labels = qc_calloc(ham->n ? ham->n : 1, sizeof(char *));
    for (size_t t = 0; t < ham->n; t++) {
        char *lab = qc_calloc((size_t)nphys + 1, 1);
        for (int p = 0; p < nphys; p++) lab[p] = 'I';
        for (int l = 0; l < lay->nl; l++) {
            int xb = (ham->t[t].x >> l) & 1, zb = (ham->t[t].z >> l) & 1;
            char c = (xb && zb) ? 'Y' : (xb ? 'X' : (zb ? 'Z' : 0));
            if (c) lab[nphys - 1 - lay->final_pos[l]] = c;
        }
        labels[t] = lab;
    }
    return labels;
}
