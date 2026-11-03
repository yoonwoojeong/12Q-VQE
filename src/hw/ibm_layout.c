/*
 * Hardware-aware layout and routing.
 *
 * layout_select: greedy growth of a connected subgraph from every seed qubit,
 *   always adding the neighbour with the smallest added error (readout error
 *   of the qubit plus CZ errors to already selected neighbours). The subgraph
 *   with the smallest total error wins.
 * layout_compile: maps the logical UCCSD circuit onto that subgraph, inserts
 *   SWAPs (3 CX each) along shortest paths inside the subgraph before every
 *   non-adjacent CX, and rewrites to the native gate set rz/sx/x/cz.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DEFAULT_READOUT 0.02
#define DEFAULT_CZ 0.01

static double qerr(const IBMCalibration *cal, int q)
{
    double r = cal->readout_err[q];
    return isnan(r) ? DEFAULT_READOUT : r;
}

static double eerr(const IBMCalibration *cal, int e)
{
    double c = cal->edge_err[e];
    return isnan(c) ? DEFAULT_CZ : c;
}

static int edge_index(const IBMCalibration *cal, int a, int b)
{
    for (int e = 0; e < cal->nedge; e++)
        if ((cal->edge_a[e] == a && cal->edge_b[e] == b) || (cal->edge_a[e] == b && cal->edge_b[e] == a)) return e;
    return -1;
}

int layout_select(const IBMCalibration *cal, int nlogical, IBMLayout *out)
{
    memset(out, 0, sizeof *out);
    out->nl = nlogical;
    if (!cal->has_cz) {
        qc_warn("backend %s has no 'cz' in its basis gates (ECR-only backends are not supported)", cal->backend_name);
        return -1;
    }
    int n = cal->nqubit;
    int *deg = qc_calloc((size_t)n, sizeof(int));
    for (int e = 0; e < cal->nedge; e++) { deg[cal->edge_a[e]]++; deg[cal->edge_b[e]]++; }

    double best = 1e300;
    int best_order[QC_MAX_QUBIT_SV], best_set[IBM_MAX_QUBIT], best_n = 0;
    int found = 0;
    for (int seed = 0; seed < n; seed++) {
        if (deg[seed] == 0) continue;
        int in[IBM_MAX_QUBIT] = {0};
        int order[QC_MAX_QUBIT_SV];
        int ns = 0;
        in[seed] = 1;
        order[ns++] = seed;
        int ok = 1;
        while (ns < nlogical) {
            int bestc = -1;
            double bestcost = 1e300;
            for (int s = 0; s < ns; s++) {
                int u = order[s];
                for (int e = 0; e < cal->nedge; e++) {
                    int c = -1;
                    if (cal->edge_a[e] == u) c = cal->edge_b[e];
                    else if (cal->edge_b[e] == u) c = cal->edge_a[e];
                    if (c < 0 || in[c]) continue;
                    /* added cost: readout of c + CZ errors to all selected neighbours of c */
                    double cost = qerr(cal, c);
                    for (int e2 = 0; e2 < cal->nedge; e2++) {
                        int d = -1;
                        if (cal->edge_a[e2] == c) d = cal->edge_b[e2];
                        else if (cal->edge_b[e2] == c) d = cal->edge_a[e2];
                        if (d >= 0 && in[d]) cost += eerr(cal, e2);
                    }
                    if (cost < bestcost || (cost == bestcost && c < bestc)) { bestcost = cost; bestc = c; }
                }
            }
            if (bestc < 0) { ok = 0; break; }
            in[bestc] = 1;
            order[ns++] = bestc;
        }
        if (!ok) continue;
        double score = 0.0;
        for (int s = 0; s < ns; s++) score += qerr(cal, order[s]);
        for (int e = 0; e < cal->nedge; e++)
            if (in[cal->edge_a[e]] && in[cal->edge_b[e]]) score += eerr(cal, e);
        if (score < best) {
            best = score;
            found = 1;
            best_n = ns;
            memcpy(best_order, order, sizeof(int) * (size_t)ns);
            memcpy(best_set, in, sizeof(int) * IBM_MAX_QUBIT);
        }
    }
    free(deg);
    if (!found) return -1;
    for (int l = 0; l < nlogical; l++) out->phys[l] = best_order[l];
    out->nsubset = 0;
    for (int q = 0; q < n; q++) if (best_set[q]) out->subset[out->nsubset++] = q;
    out->score = best;
    out->ok = 1;
    (void)best_n;
    return 0;
}

/* BFS shortest path inside the layout subgraph from a to b; returns length in path[] */
static int bfs_path(const IBMCalibration *cal, const IBMLayout *lay, int a, int b, int *path)
{
    int n = cal->nqubit;
    int *prev = qc_calloc((size_t)n, sizeof(int));
    int *queue = qc_calloc((size_t)n, sizeof(int));
    for (int i = 0; i < n; i++) prev[i] = -2;
    int in_sub[IBM_MAX_QUBIT] = {0};
    for (int i = 0; i < lay->nsubset; i++) in_sub[lay->subset[i]] = 1;
    int qh = 0, qt = 0;
    queue[qt++] = a;
    prev[a] = -1;
    while (qh < qt) {
        int u = queue[qh++];
        if (u == b) break;
        for (int e = 0; e < cal->nedge; e++) {
            int v = -1;
            if (cal->edge_a[e] == u) v = cal->edge_b[e];
            else if (cal->edge_b[e] == u) v = cal->edge_a[e];
            if (v < 0 || !in_sub[v] || prev[v] != -2) continue;
            prev[v] = u;
            queue[qt++] = v;
        }
    }
    int len = 0;
    if (prev[b] != -2) {
        int tmp[IBM_MAX_QUBIT];
        int k = 0;
        for (int v = b; v != -1; v = prev[v]) tmp[k++] = v;
        for (int i = 0; i < k; i++) path[i] = tmp[k - 1 - i];
        len = k;
    }
    free(prev);
    free(queue);
    return len;
}

int layout_compile(const Circuit *logical, const IBMCalibration *cal, IBMLayout *lay, Circuit *native_phys)
{
    int n = cal->nqubit;
    int pos[QC_MAX_QUBIT_SV];
    int at[IBM_MAX_QUBIT];
    for (int p = 0; p < n; p++) at[p] = -1;
    for (int l = 0; l < lay->nl; l++) {
        pos[l] = lay->phys[l];
        at[lay->phys[l]] = l;
    }
    Circuit phys;
    circ_init(&phys, n);
    lay->nswaps = 0;
    int path[IBM_MAX_QUBIT];
    for (int k = 0; k < logical->ng; k++) {
        const Gate *g = &logical->g[k];
        if (g->kind != G_CX) {
            if (g->kind == G_CZ) qc_die("layout_compile: logical circuit must not contain CZ");
            circ_add(&phys, g->kind, pos[g->q0], -1, g->theta);
            continue;
        }
        int c = g->q0, t = g->q1;
        for (;;) {
            int e = edge_index(cal, pos[c], pos[t]);
            if (e >= 0) break;
            int len = bfs_path(cal, lay, pos[c], pos[t], path);
            if (len < 2) qc_die("layout_compile: no route between physical %d and %d inside the layout", pos[c], pos[t]);
            int p1 = path[0], p2 = path[1];
            /* SWAP(p1,p2) = CX(p1,p2) CX(p2,p1) CX(p1,p2) */
            circ_add(&phys, G_CX, p1, p2, 0);
            circ_add(&phys, G_CX, p2, p1, 0);
            circ_add(&phys, G_CX, p1, p2, 0);
            lay->nswaps++;
            int lc = at[p1], lo = at[p2];
            at[p1] = lo; at[p2] = lc;
            if (lc >= 0) pos[lc] = p2;
            if (lo >= 0) pos[lo] = p1;
        }
        circ_add(&phys, G_CX, pos[c], pos[t], 0);
    }
    for (int l = 0; l < lay->nl; l++) lay->final_pos[l] = pos[l];

    circ_to_native(&phys, native_phys);
    circ_free(&phys);

    /* every native CZ must sit on a coupler */
    int bad = 0;
    for (int k = 0; k < native_phys->ng; k++) {
        const Gate *g = &native_phys->g[k];
        if (g->kind == G_CZ && edge_index(cal, g->q0, g->q1) < 0) bad++;
    }
    if (bad) qc_die("layout_compile: %d CZ gates are not on couplers (internal error)", bad);
    lay->ok = 1;
    return 0;
}
