/*
 * Controlled experiment: configuration file and construction of the arms
 * (bare_hf, matched_hf, vqe). See experiment.c for the commands.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void exp_config_default(ExpConfig *c)
{
    memset(c, 0, sizeof *c);
    snprintf(c->name, sizeof c->name, "h2_cas26");
    c->bond = 0.735;
    c->nel = 2;
    c->norb = 6;
    snprintf(c->backend, sizeof c->backend, "auto");
    c->ansatz_uccsd = 0;
    c->gtrunc = 1e-4;
    c->precision = 0.02;
    c->resilience = 1;
    c->seed = 12345;
    c->arms = EXP_ARM_BARE | EXP_ARM_MATCHED | EXP_ARM_VQE;
}

static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1])) *--e = 0;
    return s;
}

/* key = value lines, '#' starts a comment. Unknown keys are errors so typos are not ignored. */
int exp_config_load(const char *path, ExpConfig *c, char *err, size_t errlen)
{
    exp_config_default(c);
    snprintf(c->path, sizeof c->path, "%s", path);
    size_t len = 0;
    char *txt = qc_read_file(path, &len);
    if (!txt) { snprintf(err, errlen, "cannot read config '%s'", path); return -1; }
    int lineno = 0;
    char *save = NULL;
    int rc = 0;
    for (char *line = strtok_r(txt, "\n", &save); line; line = strtok_r(NULL, "\n", &save)) {
        lineno++;
        char *hash = strchr(line, '#');
        if (hash) *hash = 0;
        char *l = trim(line);
        if (!*l) continue;
        char *eq = strchr(l, '=');
        if (!eq) { snprintf(err, errlen, "%s:%d: expected 'key = value'", path, lineno); rc = -1; break; }
        *eq = 0;
        char *k = trim(l), *v = trim(eq + 1);
        if (strcmp(k, "name") == 0) snprintf(c->name, sizeof c->name, "%s", v);
        else if (strcmp(k, "bond_angstrom") == 0) c->bond = atof(v);
        else if (strcmp(k, "active_electrons") == 0) c->nel = atoi(v);
        else if (strcmp(k, "active_orbitals") == 0) c->norb = atoi(v);
        else if (strcmp(k, "backend") == 0) snprintf(c->backend, sizeof c->backend, "%s", v);
        else if (strcmp(k, "ansatz") == 0) {
            if (strcmp(v, "givens") == 0) c->ansatz_uccsd = 0;
            else if (strcmp(v, "uccsd") == 0) c->ansatz_uccsd = 1;
            else { snprintf(err, errlen, "%s:%d: ansatz must be givens or uccsd", path, lineno); rc = -1; break; }
        } else if (strcmp(k, "givens_trunc") == 0) c->gtrunc = atof(v);
        else if (strcmp(k, "precision") == 0) c->precision = atof(v);
        else if (strcmp(k, "resilience_level") == 0) c->resilience = atoi(v);
        else if (strcmp(k, "seed") == 0) c->seed = atoi(v);
        else if (strcmp(k, "arms") == 0) {
            c->arms = 0;
            char buf[256];
            snprintf(buf, sizeof buf, "%s", v);
            char *s2 = NULL;
            for (char *a = strtok_r(buf, ", ", &s2); a; a = strtok_r(NULL, ", ", &s2)) {
                if (strcmp(a, "bare_hf") == 0) c->arms |= EXP_ARM_BARE;
                else if (strcmp(a, "matched_hf") == 0) c->arms |= EXP_ARM_MATCHED;
                else if (strcmp(a, "vqe") == 0) c->arms |= EXP_ARM_VQE;
                else { snprintf(err, errlen, "%s:%d: unknown arm '%s'", path, lineno, a); rc = -1; break; }
            }
            if (rc) break;
        } else { snprintf(err, errlen, "%s:%d: unknown key '%s'", path, lineno, k); rc = -1; break; }
    }
    free(txt);
    if (rc == 0) {
        if (!(c->arms & EXP_ARM_VQE)) { snprintf(err, errlen, "arms must include vqe"); rc = -1; }
        else if (c->precision <= 0) { snprintf(err, errlen, "precision must be > 0"); rc = -1; }
        else if (c->resilience < 0 || c->resilience > 2) { snprintf(err, errlen, "resilience_level must be 0, 1 or 2"); rc = -1; }
        else if (c->nel != 2 || 2 * c->norb > QC_MAX_QUBIT_SV) {
            snprintf(err, errlen, "this build supports 2 active electrons in at most %d orbitals", QC_MAX_QUBIT_SV / 2);
            rc = -1;
        }
    }
    return rc;
}

/* Logical circuits of the enabled arms (order: bare_hf, matched_hf, vqe) and their exact energies. */
int exp_build_arms(const ExpConfig *cfg, const UCCSD *u, const double *theta, const PauliSum *ham,
                   ExpArm *arms, int *narms, double *vqe_fidelity)
{
    StateVec target, v;
    sv_init(&target, u->nq);
    sv_init(&v, u->nq);
    uccsd_prepare(u, theta, &target);
    Circuit treat;
    *vqe_fidelity = 1.0;
    if (cfg->ansatz_uccsd) {
        uccsd_circuit(u, theta, &treat);
    } else {
        GivensInfo gi;
        givens_state_circuit(&target, u->nq, cfg->gtrunc, &treat, &gi);
        *vqe_fidelity = gi.fidelity;
        if (gi.fidelity < 0.999) qc_die("Givens circuit fidelity %.6f is below 0.999; lower givens_trunc", gi.fidelity);
    }
    int n = 0;
    if (cfg->arms & EXP_ARM_BARE) {
        ExpArm *a = &arms[n++];
        a->name = "bare_hf";
        a->purpose = "X gates only: readout control";
        circ_init(&a->logical, u->nq);
        for (int q = 0; q < u->nq; q++)
            if ((u->hf_mask >> q) & 1) circ_add(&a->logical, G_X, q, -1, 0);
    }
    if (cfg->arms & EXP_ARM_MATCHED) {
        ExpArm *a = &arms[n++];
        a->name = "matched_hf";
        a->purpose = "same gates, angles 0: noise control";
        circ_init(&a->logical, u->nq);
        for (int k = 0; k < treat.ng; k++) {
            Gate g = treat.g[k];
            if (g.kind == G_RZ) g.theta = 0.0;
            circ_add(&a->logical, g.kind, g.q0, g.q1, g.theta);
        }
    }
    {
        ExpArm *a = &arms[n++];
        a->name = "vqe";
        a->purpose = cfg->ansatz_uccsd ? "UCCSD circuit: treatment" : "Givens VQE state: treatment";
        a->logical = treat;
    }
    for (int i = 0; i < n; i++) {
        sv_set_basis(&v, 0);
        circ_simulate(&arms[i].logical, &v);
        arms[i].exact = sv_expect_sum(&v, ham);
    }
    *narms = n;
    sv_free(&target);
    sv_free(&v);
    return 0;
}

