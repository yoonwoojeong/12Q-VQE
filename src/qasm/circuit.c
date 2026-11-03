/*
 * Gate-level circuit IR, Pauli-rotation compilation, native decomposition,
 * simulation of the IR, and OpenQASM 3 emission.
 *
 * Pauli rotation exp(i phi P), P = sigma_{q1} x ... (Hermitian):
 *   basis change: X -> H ; Y -> SDG then H  (maps the Pauli to Z)
 *   nearest-neighbour CNOT chain along the support, RZ(-2 phi) on the last qubit, undo.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <complex.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

void circ_init(Circuit *c, int nq)
{
    c->nq = nq;
    c->ng = 0;
    c->cap = 256;
    c->g = qc_calloc((size_t)c->cap, sizeof(Gate));
}

void circ_free(Circuit *c)
{
    free(c->g);
    memset(c, 0, sizeof *c);
}

void circ_add(Circuit *c, GateKind k, int q0, int q1, double theta)
{
    if (c->ng == c->cap) {
        c->cap *= 2;
        c->g = qc_realloc(c->g, (size_t)c->cap * sizeof(Gate));
    }
    Gate g = {k, q0, q1, theta};
    c->g[c->ng++] = g;
}

void circ_append_pauli_rotation(Circuit *c, uint32_t x, uint32_t z, double phi, int nq)
{
    (void)nq;
    uint32_t supp = x | z;
    if (supp == 0) return;   /* identity: global phase only */
    int qs[QC_MAX_QUBIT_SV];
    int nsup = 0;
    for (int q = 0; q < 32; q++)
        if ((supp >> q) & 1) qs[nsup++] = q;
    int t = qs[nsup - 1];
    (void)t;
    /* basis change into Z */
    for (int k = 0; k < nsup; k++) {
        int q = qs[k];
        int xb = (x >> q) & 1, zb = (z >> q) & 1;
        if (xb && !zb) circ_add(c, G_H, q, -1, 0);
        else if (xb && zb) { circ_add(c, G_SDG, q, -1, 0); circ_add(c, G_H, q, -1, 0); }
    }
    /* nearest-neighbour parity chain: CX(q0,q1) CX(q1,q2) ... leaves the parity on the last
       qubit and only couples consecutive support qubits (JW strings are contiguous) */
    for (int k = 0; k < nsup - 1; k++) circ_add(c, G_CX, qs[k], qs[k + 1], 0);
    circ_add(c, G_RZ, qs[nsup - 1], -1, -2.0 * phi);   /* exp(i phi Z) = RZ(-2 phi) */
    for (int k = nsup - 2; k >= 0; k--) circ_add(c, G_CX, qs[k], qs[k + 1], 0);
    for (int k = nsup - 1; k >= 0; k--) {
        int q = qs[k];
        int xb = (x >> q) & 1, zb = (z >> q) & 1;
        if (xb && !zb) circ_add(c, G_H, q, -1, 0);
        else if (xb && zb) { circ_add(c, G_H, q, -1, 0); circ_add(c, G_S, q, -1, 0); }
    }
}

/* ---- native decomposition: rz, sx, x, cz ---- */

static void emit_h_native(Circuit *o, int q)
{
    /* H = e^{i pi/2} RZ(pi/2) SX RZ(pi/2) (global phase dropped) */
    circ_add(o, G_RZ, q, -1, QC_PI / 2);
    circ_add(o, G_SX, q, -1, 0);
    circ_add(o, G_RZ, q, -1, QC_PI / 2);
}

void circ_to_native(const Circuit *in, Circuit *out)
{
    circ_init(out, in->nq);
    for (int k = 0; k < in->ng; k++) {
        const Gate *g = &in->g[k];
        switch (g->kind) {
        case G_H: emit_h_native(out, g->q0); break;
        case G_SDG: circ_add(out, G_RZ, g->q0, -1, -QC_PI / 2); break;
        case G_S: circ_add(out, G_RZ, g->q0, -1, QC_PI / 2); break;
        case G_X: circ_add(out, G_X, g->q0, -1, 0); break;
        case G_SX: circ_add(out, G_SX, g->q0, -1, 0); break;
        case G_RZ: circ_add(out, G_RZ, g->q0, -1, g->theta); break;
        case G_CZ: circ_add(out, G_CZ, g->q0, g->q1, 0); break;
        case G_CX:
            emit_h_native(out, g->q1);
            circ_add(out, G_CZ, g->q0, g->q1, 0);
            emit_h_native(out, g->q1);
            break;
        default: qc_die("circ_to_native: unsupported gate kind %d", (int)g->kind);
        }
    }
}

/* ---- statistics ---- */

int circ_count_2q(const Circuit *c)
{
    int n = 0;
    for (int k = 0; k < c->ng; k++)
        if (c->g[k].kind == G_CX || c->g[k].kind == G_CZ) n++;
    return n;
}

int circ_count_kind(const Circuit *c, GateKind kind)
{
    int n = 0;
    for (int k = 0; k < c->ng; k++)
        if (c->g[k].kind == kind) n++;
    return n;
}

int circ_depth(const Circuit *c)
{
    int *lay = qc_calloc((size_t)c->nq, sizeof(int));
    int depth = 0;
    for (int k = 0; k < c->ng; k++) {
        const Gate *g = &c->g[k];
        int a = g->q0, b = g->q1;
        int lvl = lay[a];
        if (b >= 0 && lay[b] > lvl) lvl = lay[b];
        lvl++;
        lay[a] = lvl;
        if (b >= 0) lay[b] = lvl;
        if (lvl > depth) depth = lvl;
    }
    free(lay);
    return depth;
}

/* ---- simulation of the IR (validation) ---- */

void circ_simulate(const Circuit *c, StateVec *sv)
{
    static const double complex s2 = 0.70710678118654752440;
    for (int k = 0; k < c->ng; k++) {
        const Gate *g = &c->g[k];
        double complex u[2][2];
        switch (g->kind) {
        case G_H: u[0][0] = s2; u[0][1] = s2; u[1][0] = s2; u[1][1] = -s2; sv_apply_1q(sv, g->q0, u); break;
        case G_X: u[0][0] = 0; u[0][1] = 1; u[1][0] = 1; u[1][1] = 0; sv_apply_1q(sv, g->q0, u); break;
        case G_S: u[0][0] = 1; u[0][1] = 0; u[1][0] = 0; u[1][1] = I; sv_apply_1q(sv, g->q0, u); break;
        case G_SDG: u[0][0] = 1; u[0][1] = 0; u[1][0] = 0; u[1][1] = -I; sv_apply_1q(sv, g->q0, u); break;
        case G_SX:
            u[0][0] = 0.5 * (1 + I); u[0][1] = 0.5 * (1 - I);
            u[1][0] = 0.5 * (1 - I); u[1][1] = 0.5 * (1 + I);
            sv_apply_1q(sv, g->q0, u);
            break;
        case G_RZ:
            u[0][0] = cexp(-I * g->theta / 2); u[0][1] = 0;
            u[1][0] = 0; u[1][1] = cexp(I * g->theta / 2);
            sv_apply_1q(sv, g->q0, u);
            break;
        case G_CX: sv_apply_cx(sv, g->q0, g->q1); break;
        case G_CZ: sv_apply_cz(sv, g->q0, g->q1); break;
        default: qc_die("circ_simulate: unsupported gate");
        }
    }
}

/* ---- OpenQASM 3 ---- */

static void sput(char **buf, size_t *len, size_t *cap, const char *s)
{
    size_t n = strlen(s);
    if (*len + n + 1 > *cap) {
        *cap = (*len + n + 1) * 2 + 256;
        *buf = qc_realloc(*buf, *cap);
    }
    memcpy(*buf + *len, s, n + 1);
    *len += n;
}

static void sputf(char **buf, size_t *len, size_t *cap, const char *fmt, ...)
{
    char tmp[256];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(tmp, sizeof tmp, fmt, ap);
    va_end(ap);
    sput(buf, len, cap, tmp);
}

char *qasm3_string(const Circuit *c, int nq_declared, int native)
{
    (void)native;   /* the IR is already native when produced by circ_to_native */
    char *buf = NULL;
    size_t len = 0, cap = 0;
    sput(&buf, &len, &cap, "OPENQASM 3.0;\ninclude \"stdgates.inc\";\n");
    sputf(&buf, &len, &cap, "qubit[%d] q;\n", nq_declared > c->nq ? nq_declared : c->nq);
    for (int k = 0; k < c->ng; k++) {
        const Gate *g = &c->g[k];
        switch (g->kind) {
        case G_H: sputf(&buf, &len, &cap, "h q[%d];\n", g->q0); break;
        case G_X: sputf(&buf, &len, &cap, "x q[%d];\n", g->q0); break;
        case G_S: sputf(&buf, &len, &cap, "s q[%d];\n", g->q0); break;
        case G_SDG: sputf(&buf, &len, &cap, "sdg q[%d];\n", g->q0); break;
        case G_SX: sputf(&buf, &len, &cap, "sx q[%d];\n", g->q0); break;
        case G_RZ: sputf(&buf, &len, &cap, "rz(%.17g) q[%d];\n", g->theta, g->q0); break;
        case G_CX: sputf(&buf, &len, &cap, "cx q[%d], q[%d];\n", g->q0, g->q1); break;
        case G_CZ: sputf(&buf, &len, &cap, "cz q[%d], q[%d];\n", g->q0, g->q1); break;
        default: break;
        }
    }
    return buf;
}

int qasm3_write(const Circuit *c, FILE *f, const char *const *phys_names, int native)
{
    (void)phys_names;
    char *s = qasm3_string(c, c->nq, native);
    fputs(s, f);
    free(s);
    return 0;
}
