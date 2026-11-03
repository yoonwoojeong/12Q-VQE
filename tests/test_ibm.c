/*
 * IBM layer tests on SYNTHETIC devices (in-memory fixtures, not hardware data):
 *  - calibration JSON parsing,
 *  - layout selection avoids noisy qubits and returns a connected subgraph,
 *  - routing + native compilation of the UCCSD circuit: every CZ on a coupler,
 *    and the routed state equals the logical state up to the final permutation.
 */
#include "testlib.h"
#include <cjson/cJSON.h>

/* build cfg/props JSON for a device with given couplers and per-qubit/edge errors */
static void make_device(int nq, int nedge, const int (*edges)[2], double readout_base, double cz_base,
                        char **cfg_out, char **props_out, int noisy_qubit)
{
    cJSON *cfg = cJSON_CreateObject();
    cJSON_AddStringToObject(cfg, "backend_name", "synthetic_device");
    cJSON_AddStringToObject(cfg, "backend_version", "0.0.0-fixture");
    cJSON_AddNumberToObject(cfg, "n_qubits", nq);
    cJSON *bg = cJSON_AddArrayToObject(cfg, "basis_gates");
    cJSON_AddItemToArray(bg, cJSON_CreateString("cz"));
    cJSON_AddItemToArray(bg, cJSON_CreateString("rz"));
    cJSON_AddItemToArray(bg, cJSON_CreateString("sx"));
    cJSON_AddItemToArray(bg, cJSON_CreateString("x"));
    cJSON *cm = cJSON_AddArrayToObject(cfg, "coupling_map");
    for (int e = 0; e < nedge; e++) {
        for (int d = 0; d < 2; d++) {
            cJSON *p = cJSON_CreateArray();
            cJSON_AddItemToArray(p, cJSON_CreateNumber(d ? edges[e][1] : edges[e][0]));
            cJSON_AddItemToArray(p, cJSON_CreateNumber(d ? edges[e][0] : edges[e][1]));
            cJSON_AddItemToArray(cm, p);
        }
    }
    cJSON *pr = cJSON_CreateObject();
    cJSON *qs = cJSON_AddArrayToObject(pr, "qubits");
    for (int q = 0; q < nq; q++) {
        cJSON *arr = cJSON_CreateArray();
        double ro = readout_base + 0.001 * (q % 5) + (q == noisy_qubit ? 0.25 : 0.0);
        cJSON *a = cJSON_CreateObject();
        cJSON_AddStringToObject(a, "name", "readout_error"); cJSON_AddNumberToObject(a, "value", ro);
        cJSON_AddItemToArray(arr, a);
        cJSON *t1 = cJSON_CreateObject();
        cJSON_AddStringToObject(t1, "name", "T1"); cJSON_AddNumberToObject(t1, "value", 150.0 + q);
        cJSON_AddItemToArray(arr, t1);
        cJSON *t2 = cJSON_CreateObject();
        cJSON_AddStringToObject(t2, "name", "T2"); cJSON_AddNumberToObject(t2, "value", 90.0 + q);
        cJSON_AddItemToArray(arr, t2);
        cJSON_AddItemToArray(qs, arr);
    }
    cJSON *gs = cJSON_AddArrayToObject(pr, "gates");
    for (int e = 0; e < nedge; e++) {
        cJSON *g = cJSON_CreateObject();
        cJSON_AddStringToObject(g, "gate", "cz");
        cJSON *gq = cJSON_AddArrayToObject(g, "qubits");
        char a[16], b[16];
        snprintf(a, sizeof a, "%d", edges[e][0]); snprintf(b, sizeof b, "%d", edges[e][1]);
        cJSON_AddItemToArray(gq, cJSON_CreateString(a));
        cJSON_AddItemToArray(gq, cJSON_CreateString(b));
        cJSON *ps = cJSON_AddArrayToObject(g, "parameters");
        cJSON *pe = cJSON_CreateObject();
        cJSON_AddStringToObject(pe, "name", "gate_error");
        cJSON_AddNumberToObject(pe, "value", cz_base + 0.0005 * (e % 4));
        cJSON_AddItemToArray(ps, pe);
        cJSON *pl = cJSON_CreateObject();
        cJSON_AddStringToObject(pl, "name", "gate_length");
        cJSON_AddNumberToObject(pl, "value", 60.0 + e % 3);
        cJSON_AddItemToArray(ps, pl);
        cJSON_AddItemToArray(gs, g);
    }
    *cfg_out = cJSON_PrintUnformatted(cfg);
    *props_out = cJSON_PrintUnformatted(pr);
    cJSON_Delete(cfg);
    cJSON_Delete(pr);
}

int main(void)
{
    printf("== test_ibm (synthetic fixtures)\n");

    /* --- ring of 12 qubits: routing + equivalence on a statevector --- */
    int ring[12][2];
    for (int i = 0; i < 12; i++) { ring[i][0] = i; ring[i][1] = (i + 1) % 12; }
    char *cfg, *props;
    make_device(12, 12, ring, 0.01, 0.004, &cfg, &props, -1);
    IBMCalibration cal;
    CHECK(ibm_calibration_parse(cfg, props, &cal) == 0, "calibration JSON parsed");
    CHECK(cal.nqubit == 12 && cal.nedge == 12 && cal.has_cz, "12 qubits, 12 couplers, cz basis");
    CHECK(fabs(cal.edge_err[0] - 0.004) < 1e-12 && fabs(cal.readout_err[3] - 0.013) < 1e-12,
          "gate and readout errors read from properties");
    IBMLayout lay;
    CHECK(layout_select(&cal, 12, &lay) == 0 && lay.nsubset == 12, "layout selected 12 physical qubits");

    Sys s;
    sys_build(&s, 2, 6);
    UCCSD u;
    uccsd_build(&u, 12, 2);
    double *th = qc_calloc((size_t)u.nexc, sizeof(double));
    unsigned seed = 11;
    for (int k = 0; k < u.nexc; k++) { seed = seed * 1103515245u + 12345u; th[k] = (((seed >> 8) & 0xffff) / 65535.0 - 0.5); }
    Circuit logical, native;
    uccsd_circuit(&u, th, &logical);
    layout_compile(&logical, &cal, &lay, &native);
    CHECK(lay.nswaps > 0, "ring topology needs SWAP routing (%d SWAPs inserted)", lay.nswaps);

    /* every CZ of the native circuit is on a coupler */
    int offcoupler = 0;
    for (int k = 0; k < native.ng; k++) {
        const Gate *g = &native.g[k];
        if (g->kind != G_CZ) continue;
        int ok = 0;
        for (int e = 0; e < cal.nedge; e++)
            if ((cal.edge_a[e] == g->q0 && cal.edge_b[e] == g->q1) || (cal.edge_a[e] == g->q1 && cal.edge_b[e] == g->q0)) ok = 1;
        if (!ok) offcoupler++;
    }
    CHECK(offcoupler == 0, "all native CZ gates act on couplers (violations %d)", offcoupler);
    int only_native = 1;
    for (int k = 0; k < native.ng; k++) {
        GateKind kd = native.g[k].kind;
        if (!(kd == G_RZ || kd == G_SX || kd == G_X || kd == G_CZ)) only_native = 0;
    }
    CHECK(only_native, "native circuit uses only rz, sx, x, cz");

    /* routed physical state == logical state, up to the final logical->physical permutation */
    StateVec L, P;
    sv_init(&L, 12);
    sv_init(&P, 12);
    circ_simulate(&logical, &L);
    circ_simulate(&native, &P);     /* physical register starts in |0...0> like the logical one */
    double maxdev = 0.0;
    double complex phase = 0.0;
    /* overlap check with the permutation: amplitude of logical i sits at physical j(i) */
    double complex ov = 0.0;
    for (size_t i = 0; i < L.dim; i++) {
        size_t j = 0;
        for (int l = 0; l < 12; l++) if ((i >> l) & 1) j |= (size_t)1 << lay.final_pos[l];
        ov += conj(L.a[i]) * P.a[j];
        (void)phase;
    }
    double ovmag = cabs(ov);
    for (size_t i = 0; i < L.dim; i++) {
        size_t j = 0;
        for (int l = 0; l < 12; l++) if ((i >> l) & 1) j |= (size_t)1 << lay.final_pos[l];
        double d = cabs(P.a[j] - L.a[i] * (ov / cabs(ov)));
        if (d > maxdev) maxdev = d;
    }
    CHECK(fabs(ovmag - 1.0) < 1e-9, "routed state equals logical UCCSD state up to permutation (|<.|.>| = %.12f)", ovmag);
    CHECK(maxdev < 1e-9, "amplitude-by-amplitude agreement after permutation (max dev %.2e)", maxdev);
    /* Givens circuit on the ring: contiguous support, so no SWAPs, and equivalent after routing */
    {
        StateVec tgt;
        sv_init(&tgt, 12);
        uccsd_prepare(&u, th, &tgt);
        Circuit gl, gn;
        GivensInfo gi;
        givens_state_circuit(&tgt, 12, 0.0, &gl, &gi);
        IBMLayout glay;
        layout_select(&cal, 12, &glay);
        layout_compile(&gl, &cal, &glay, &gn);
        CHECK(glay.nswaps == 0, "Givens circuit routes on a path layout with no SWAPs (%d)", glay.nswaps);
        StateVec P2;
        sv_init(&P2, 12);
        circ_simulate(&gn, &P2);
        double complex ov2 = 0.0;
        for (size_t i = 0; i < tgt.dim; i++) {
            size_t jj = 0;
            for (int l = 0; l < 12; l++) if ((i >> l) & 1) jj |= (size_t)1 << glay.final_pos[l];
            ov2 += conj(tgt.a[i]) * P2.a[jj];
        }
        CHECK(fabs(cabs(ov2) - 1.0) < 1e-9, "routed Givens circuit prepares the UCCSD state (|<.|.>| = %.12f)", cabs(ov2));
        double dur = circ_duration_ns(&gn, &cal);
        CHECK(dur > 0, "timed critical path computed (%.1f us)", dur / 1000.0);
        circ_free(&gl); circ_free(&gn); sv_free(&tgt); sv_free(&P2);
    }
    free(cfg); free(props);

    /* --- 4x4 grid with a noisy corner: layout must avoid the noisy qubit --- */
    int grid_e[24][2];
    int ne = 0;
    for (int r = 0; r < 4; r++)
        for (int c = 0; c < 4; c++) {
            int q = r * 4 + c;
            if (c < 3) { grid_e[ne][0] = q; grid_e[ne][1] = q + 1; ne++; }
            if (r < 3) { grid_e[ne][0] = q; grid_e[ne][1] = q + 4; ne++; }
        }
    make_device(16, ne, grid_e, 0.01, 0.005, &cfg, &props, 0);   /* qubit 0 is very noisy */
    IBMCalibration gcal;
    CHECK(ibm_calibration_parse(cfg, props, &gcal) == 0 && gcal.nedge == ne, "grid device parsed (%d couplers)", gcal.nedge);
    IBMLayout gl;
    CHECK(layout_select(&gcal, 12, &gl) == 0, "layout found on 4x4 grid");
    int has_noisy = 0;
    for (int i = 0; i < gl.nsubset; i++) if (gl.subset[i] == 0) has_noisy = 1;
    CHECK(!has_noisy, "noisy qubit 0 (readout +0.25) excluded from the layout");
    Circuit gnat;
    layout_compile(&logical, &gcal, &gl, &gnat);
    int gbad = 0;
    for (int k = 0; k < gnat.ng; k++) {
        const Gate *g = &gnat.g[k];
        if (g->kind != G_CZ) continue;
        int ok = 0;
        for (int e = 0; e < gcal.nedge; e++)
            if ((gcal.edge_a[e] == g->q0 && gcal.edge_b[e] == g->q1) || (gcal.edge_a[e] == g->q1 && gcal.edge_b[e] == g->q0)) ok = 1;
        if (!ok) gbad++;
    }
    CHECK(gbad == 0, "grid-compiled CZ gates all on couplers (%d SWAPs, %d native CZ)", gl.nswaps, circ_count_kind(&gnat, G_CZ));

    free(cfg); free(props);
    free(th);
    circ_free(&logical); circ_free(&native); circ_free(&gnat);
    uccsd_free(&u);
    sv_free(&L); sv_free(&P);
    fermion_free(&s.fh);
    active_space_free(&s.as);
    rhf_free(&s.rhf);
    ao_integrals_free(&s.ints);
    basis_free(&s.bs);

    /* base64 decoding of little-endian f64 */
    double *dv = NULL; size_t dn = 0;
    /* 1.0 and 2.5 as little-endian f64 */
    const char *b64b = "AAAAAAAA8D8AAAAAAAAEQA==";
    int rc = base64_decode_f64(b64b, &dv, &dn);
    CHECK(rc == 0 && dn == 2 && fabs(dv[0] - 1.0) < 1e-15 && fabs(dv[1] - 2.5) < 1e-15, "base64 f64 decoding (1.0, 2.5)");
    free(dv);

    printf("%s\n", g_fail ? "test_ibm: FAILED" : "test_ibm: passed");
    return g_fail ? 1 : 0;
}
