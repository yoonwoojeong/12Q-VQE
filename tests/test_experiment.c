/*
 * Controlled-experiment setup: config parsing and the construction of the arms.
 */
#define _POSIX_C_SOURCE 200809L
#include "testlib.h"
#include <unistd.h>

static int write_tmp(const char *path, const char *text)
{
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fputs(text, f);
    fclose(f);
    return 0;
}

int main(void)
{
    printf("== test_experiment\n");
    char err[256];
    ExpConfig c;

    /* the shipped config */
    CHECK(exp_config_load(QC_DATA_DIR "/../experiments/h2_cas26.ini", &c, err, sizeof err) == 0,
          "shipped config parses (%s)", c.name);
    CHECK(c.arms == (EXP_ARM_BARE | EXP_ARM_MATCHED | EXP_ARM_VQE) && c.nel == 2 && c.norb == 6 &&
          c.ansatz_uccsd == 0 && fabs(c.precision - 0.02) < 1e-15 && c.seed == 12345, "shipped config values");

    write_tmp("/tmp/qchem_exp_bad1.ini", "name = x\nprecison = 0.02\n");
    CHECK(exp_config_load("/tmp/qchem_exp_bad1.ini", &c, err, sizeof err) != 0 && strstr(err, "precison"),
          "unknown key rejected (%s)", err);
    write_tmp("/tmp/qchem_exp_bad2.ini", "arms = bare_hf, vqe, magic\n");
    CHECK(exp_config_load("/tmp/qchem_exp_bad2.ini", &c, err, sizeof err) != 0 && strstr(err, "magic"),
          "unknown arm rejected (%s)", err);
    write_tmp("/tmp/qchem_exp_bad3.ini", "arms = bare_hf\n");
    CHECK(exp_config_load("/tmp/qchem_exp_bad3.ini", &c, err, sizeof err) != 0, "config without the vqe arm rejected");
    unlink("/tmp/qchem_exp_bad1.ini");
    unlink("/tmp/qchem_exp_bad2.ini");
    unlink("/tmp/qchem_exp_bad3.ini");

    /* arms */
    Sys s;
    sys_build(&s, 2, 6);
    UCCSD u;
    uccsd_build(&u, 12, 2);
    double *th = qc_calloc((size_t)u.nexc, sizeof(double));
    for (int k = 0; k < u.nexc; k++) th[k] = 0.02 * ((k % 5) - 2);
    exp_config_default(&c);
    ExpArm arms[3];
    int n = 0;
    double fid = 0;
    exp_build_arms(&c, &u, th, &s.fh.ham, arms, &n, &fid);
    CHECK(n == 3 && strcmp(arms[0].name, "bare_hf") == 0 && strcmp(arms[1].name, "matched_hf") == 0 &&
          strcmp(arms[2].name, "vqe") == 0, "three arms in order bare_hf, matched_hf, vqe");
    double ehf = active_hf_energy(&s.as);
    CHECK(fabs(arms[0].exact - ehf) < 1e-10, "bare_hf prepares HF exactly (%.12f)", arms[0].exact);
    CHECK(fabs(arms[1].exact - ehf) < 1e-10, "matched_hf prepares HF exactly (%.12f)", arms[1].exact);
    StateVec ref;
    sv_init(&ref, 12);
    double eu = vqe_energy(&u, &s.fh.ham, th, &ref);
    sv_free(&ref);
    CHECK(fabs(arms[2].exact - eu) < 1e-3 && fabs(arms[2].exact - ehf) > 1e-3,
          "vqe arm reproduces the UCCSD state's energy (%.8f vs %.8f)", arms[2].exact, eu);
    CHECK(fid > 0.999, "treatment fidelity %.8f", fid);
    int same = arms[1].logical.ng == arms[2].logical.ng;
    for (int k = 0; same && k < arms[1].logical.ng; k++) {
        const Gate *a = &arms[1].logical.g[k], *b = &arms[2].logical.g[k];
        same = a->kind == b->kind && a->q0 == b->q0 && a->q1 == b->q1 && (a->kind != G_RZ ? a->theta == b->theta : a->theta == 0.0);
    }
    CHECK(same, "matched_hf has the gates of vqe exactly, with RZ angles 0 (%d gates)", arms[1].logical.ng);
    int only_x = 1;
    for (int k = 0; k < arms[0].logical.ng; k++) if (arms[0].logical.g[k].kind != G_X) only_x = 0;
    CHECK(only_x && arms[0].logical.ng == 2, "bare_hf is two X gates");

    /* same physical qubits for every arm on a synthetic 12-qubit line */
    int line[11][2];
    for (int i = 0; i < 11; i++) { line[i][0] = i; line[i][1] = i + 1; }
    char cfg[4096], props[8192];
    int off = snprintf(cfg, sizeof cfg, "{\"backend_name\":\"synthetic_line\",\"backend_version\":\"0.0.0\",\"n_qubits\":12,"
                                       "\"basis_gates\":[\"cz\",\"rz\",\"sx\",\"x\"],\"coupling_map\":[");
    for (int i = 0; i < 11; i++) off += snprintf(cfg + off, sizeof cfg - off, "%s[%d,%d]", i ? "," : "", line[i][0], line[i][1]);
    snprintf(cfg + off, sizeof cfg - off, "]}");
    snprintf(props, sizeof props, "{\"qubits\":[],\"gates\":[]}");
    IBMCalibration cal;
    CHECK(ibm_calibration_parse(cfg, props, &cal) == 0 && cal.nedge == 11, "synthetic 12-qubit line parsed");
    IBMLayout lay;
    layout_select(&cal, 12, &lay);
    int nsw = 0, pos_same = 1;
    IBMLayout l[3];
    Circuit nat[3];
    for (int i = 0; i < 3; i++) {
        l[i] = lay;
        layout_compile(&arms[i].logical, &cal, &l[i], &nat[i]);
        nsw += l[i].nswaps;
        for (int q = 0; q < 12; q++) if (l[i].phys[q] != lay.phys[q]) pos_same = 0;
    }
    CHECK(pos_same && nsw == 0, "all arms on the same physical qubits, no SWAPs");
    CHECK(circ_count_kind(&nat[1], G_CZ) == circ_count_kind(&nat[2], G_CZ) && circ_count_kind(&nat[0], G_CZ) == 0,
          "control and treatment have the same CZ count (%d); bare_hf has none", circ_count_kind(&nat[2], G_CZ));

    for (int i = 0; i < 3; i++) { circ_free(&arms[i].logical); circ_free(&nat[i]); }
    free(th);
    uccsd_free(&u);
    fermion_free(&s.fh);
    active_space_free(&s.as);
    rhf_free(&s.rhf);
    ao_integrals_free(&s.ints);
    basis_free(&s.bs);
    printf("%s\n", g_fail ? "test_experiment: FAILED" : "test_experiment: passed");
    return g_fail ? 1 : 0;
}
