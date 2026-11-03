/*
 * IBM Quantum Compute Service REST API (libcurl + cJSON).
 *
 * Verified against the OpenAPI spec for IBM-API-Version 2026-04-15:
 *   POST https://iam.cloud.ibm.com/identity/token            (API key -> IAM bearer token)
 *   GET  /api/v1/backends
 *   GET  /api/v1/backends/{id}/configuration|properties|status
 *   POST /api/v1/jobs              {program_id, backend, params:{pubs, version, resilience_level, options}}
 *   GET  /api/v1/jobs/{id}, /api/v1/jobs/{id}/results
 * Credentials come from IBM_QUANTUM_API_KEY and IBM_QUANTUM_CRN; they are never printed.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cjson/cJSON.h>
#include <curl/curl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define IBM_IAM_URL "https://iam.cloud.ibm.com/identity/token"
#define IBM_DEFAULT_URL "https://quantum.cloud.ibm.com/api/v1"
#define IBM_EU_DE_URL "https://eu-de.quantum.cloud.ibm.com/api/v1"
/* IBM's front end answers requests without a User-Agent with HTTP 403 and an HTML page,
   and libcurl sends none by default. */
#define IBM_USER_AGENT "qchem-vqe/" QC_VERSION " (libcurl)"
#define IBM_DEFAULT_VERSION "2026-04-15"

typedef struct {
    char *data;
    size_t len;
} Buf;

static size_t write_cb(char *ptr, size_t size, size_t nmemb, void *userdata)
{
    Buf *b = userdata;
    size_t n = size * nmemb;
    b->data = qc_realloc(b->data, b->len + n + 1);
    memcpy(b->data + b->len, ptr, n);
    b->len += n;
    b->data[b->len] = 0;
    return n;
}

/* Region of an IBM Cloud CRN: crn:v1:bluemix:public:quantum-computing:<region>:a/<account>:<instance>:: */
int ibm_crn_region(const char *crn, char *out, size_t n)
{
    out[0] = 0;
    const char *p = crn;
    for (int field = 0; field < 5; field++) {
        p = strchr(p, ':');
        if (!p) return -1;
        p++;
    }
    const char *e = strchr(p, ':');
    size_t len = e ? (size_t)(e - p) : strlen(p);
    if (len == 0 || len >= n) return -1;
    memcpy(out, p, len);
    out[len] = 0;
    return 0;
}

/* API base URL for the instance's region (only eu-de differs from the default). */
const char *ibm_url_for_crn(const char *crn)
{
    char region[64];
    if (ibm_crn_region(crn, region, sizeof region) == 0 && strcmp(region, "eu-de") == 0) return IBM_EU_DE_URL;
    return IBM_DEFAULT_URL;
}

/* Explains a failed REST call, including IBM's own error message (it holds no secrets). */
void ibm_report_http_error(const IBMSession *s, const char *method, const char *path, long code, const char *resp)
{
    fprintf(stderr, "error: %s %s returned HTTP %ld\n", method, path, code);
    if (resp && *resp) {
        char msg[601];
        size_t k = 0;
        for (const char *c = resp; *c && k < sizeof msg - 1; c++) msg[k++] = (*c == '\n' || *c == '\r') ? ' ' : *c;
        msg[k] = 0;
        fprintf(stderr, "  server said: %s%s\n", msg, strlen(resp) > k ? " ..." : "");
    }
    if (code == 401 || code == 403) {
        char region[64] = "?";
        ibm_crn_region(s->crn, region, sizeof region);
        fprintf(stderr, "  endpoint: %s (instance CRN region: %s)\n", s->base_url, region);
        fprintf(stderr, "  hints: the CRN must be a Quantum instance this API key can access (check the\n"
                        "         instance's access policies); eu-de instances use %s;\n"
                        "         re-enter credentials with `qchem-vqe auth login`.\n", IBM_EU_DE_URL);
    }
}

int ibm_session_env(IBMSession *s)
{
    memset(s, 0, sizeof *s);
    /* stored credentials first (file written by `auth login`), then environment overrides */
    char path[1024];
    if (cred_default_path(path, sizeof path) == 0) cred_load(path, s);
    qc_env_get("IBM_QUANTUM_API_KEY", s->api_key, sizeof s->api_key);
    qc_env_get("IBM_QUANTUM_CRN", s->crn, sizeof s->crn);
    qc_env_get("IBM_QUANTUM_API_URL", s->base_url, sizeof s->base_url);
    qc_env_get("IBM_QUANTUM_API_VERSION", s->api_version, sizeof s->api_version);
    if (!s->base_url[0]) snprintf(s->base_url, sizeof s->base_url, "%s", ibm_url_for_crn(s->crn));
    if (!s->api_version[0]) snprintf(s->api_version, sizeof s->api_version, "%s", IBM_DEFAULT_VERSION);
    int ok = 1;
    if (!s->api_key[0]) {
        fprintf(stderr, "error: no IBM API key (run `qchem-vqe auth login` or set IBM_QUANTUM_API_KEY)\n");
        ok = 0;
    }
    if (!s->crn[0]) {
        fprintf(stderr, "error: no IBM instance CRN (run `qchem-vqe auth login` or set IBM_QUANTUM_CRN)\n");
        ok = 0;
    }
    return ok ? 0 : -1;
}

int ibm_session_authenticate(IBMSession *s)
{
    CURL *curl = curl_easy_init();
    if (!curl) return -1;
    char *key_esc = curl_easy_escape(curl, s->api_key, 0);
    char body[1024];
    snprintf(body, sizeof body, "grant_type=urn:ibm:params:oauth:grant-type:apikey&apikey=%s", key_esc ? key_esc : "");
    curl_free(key_esc);
    Buf b = {0};
    struct curl_slist *hdr = curl_slist_append(NULL, "Content-Type: application/x-www-form-urlencoded");
    hdr = curl_slist_append(hdr, "Accept: application/json");
    curl_easy_setopt(curl, CURLOPT_URL, IBM_IAM_URL);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body);
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, hdr);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, IBM_USER_AGENT);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_cb);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &b);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 60L);
    CURLcode rc = curl_easy_perform(curl);
    long code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
    curl_slist_free_all(hdr);
    curl_easy_cleanup(curl);
    memset(body, 0, sizeof body);   /* do not keep the key around */
    if (rc != CURLE_OK || code != 200 || !b.data) {
        fprintf(stderr, "error: IAM token request failed (HTTP %ld, curl %d)\n", code, (int)rc);
        free(b.data);
        return -1;
    }
    cJSON *j = cJSON_Parse(b.data);
    free(b.data);
    if (!j) return -1;
    cJSON *tok = cJSON_GetObjectItem(j, "access_token");
    cJSON *exp = cJSON_GetObjectItem(j, "expires_in");
    if (!cJSON_IsString(tok)) { cJSON_Delete(j); return -1; }
    snprintf(s->token, sizeof s->token, "%s", tok->valuestring);
    s->token_expiry = (double)time(NULL) + (cJSON_IsNumber(exp) ? exp->valuedouble : 3600.0) - 60.0;
    cJSON_Delete(j);
    return 0;
}

int ibm_http(IBMSession *s, const char *method, const char *path, const char *body, char **resp, long *status)
{
    if (!s->token[0] || (double)time(NULL) > s->token_expiry) {
        if (ibm_session_authenticate(s) != 0) return -1;
    }
    CURL *curl = curl_easy_init();
    if (!curl) return -1;
    char url[1024];
    snprintf(url, sizeof url, "%s%s", s->base_url, path);
    char auth[4200], crn[600], ver[64];
    snprintf(auth, sizeof auth, "Authorization: Bearer %s", s->token);
    snprintf(crn, sizeof crn, "Service-CRN: %s", s->crn);
    snprintf(ver, sizeof ver, "IBM-API-Version: %s", s->api_version);
    struct curl_slist *hdr = curl_slist_append(NULL, auth);
    hdr = curl_slist_append(hdr, crn);
    hdr = curl_slist_append(hdr, ver);
    hdr = curl_slist_append(hdr, "Accept: application/json");
    if (body) hdr = curl_slist_append(hdr, "Content-Type: application/json");
    Buf b = {0};
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_CUSTOMREQUEST, method);
    if (body) curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body);
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, hdr);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, IBM_USER_AGENT);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_cb);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &b);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 120L);
    CURLcode rc = curl_easy_perform(curl);
    long code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &code);
    curl_slist_free_all(hdr);
    curl_easy_cleanup(curl);
    if (status) *status = code;
    if (resp) *resp = b.data ? b.data : qc_strdup("");
    else free(b.data);
    if (rc != CURLE_OK) {
        fprintf(stderr, "error: HTTP %s %s failed (curl %d)\n", method, path, (int)rc);
        return -1;
    }
    return 0;
}

/* ---- backends ---- */

int ibm_list_backends(IBMSession *s, IBMBackendInfo **out, int *nout)
{
    char *resp = NULL;
    long code = 0;
    *out = NULL;
    *nout = 0;
    if (ibm_http(s, "GET", "/backends", NULL, &resp, &code) != 0 || code != 200) {
        ibm_report_http_error(s, "GET", "/backends", code, resp);
        free(resp);
        return -1;
    }
    cJSON *j = cJSON_Parse(resp);
    free(resp);
    cJSON *devs = cJSON_GetObjectItem(j, "devices");
    int n = cJSON_IsArray(devs) ? cJSON_GetArraySize(devs) : 0;
    IBMBackendInfo *arr = qc_calloc((size_t)(n ? n : 1), sizeof(IBMBackendInfo));
    int k = 0;
    for (cJSON *d = devs ? devs->child : NULL; d; d = d->next) {
        IBMBackendInfo *bi = &arr[k++];
        cJSON *nm = cJSON_GetObjectItem(d, "name");
        cJSON *q = cJSON_GetObjectItem(d, "qubits");
        cJSON *sim = cJSON_GetObjectItem(d, "is_simulator");
        cJSON *st = cJSON_GetObjectItem(d, "status");
        cJSON *ql = cJSON_GetObjectItem(d, "queue_length");
        snprintf(bi->name, sizeof bi->name, "%s", cJSON_IsString(nm) ? nm->valuestring : "");
        bi->n_qubits = cJSON_IsNumber(q) ? q->valueint : 0;
        bi->is_simulator = cJSON_IsTrue(sim);
        cJSON *sn = st ? cJSON_GetObjectItem(st, "name") : NULL;
        snprintf(bi->status, sizeof bi->status, "%s", cJSON_IsString(sn) ? sn->valuestring : "unknown");
        bi->queue_length = cJSON_IsNumber(ql) ? ql->valueint : -1;
    }
    cJSON_Delete(j);
    *out = arr;
    *nout = k;
    return 0;
}

static int get_checked(IBMSession *s, const char *path, char **out)
{
    long code = 0;
    char *resp = NULL;
    if (ibm_http(s, "GET", path, NULL, &resp, &code) != 0 || code != 200) {
        ibm_report_http_error(s, "GET", path, code, resp);
        free(resp);
        if (out) *out = NULL;
        return -1;
    }
    *out = resp;
    return 0;
}

int ibm_get_backend_raw(IBMSession *s, const char *name, char **cfg_json, char **props_json, char **status_json)
{
    char path[512];
    snprintf(path, sizeof path, "/backends/%s/configuration", name);
    if (get_checked(s, path, cfg_json) != 0) return -1;
    snprintf(path, sizeof path, "/backends/%s/properties", name);
    if (get_checked(s, path, props_json) != 0) return -1;
    if (status_json) {
        snprintf(path, sizeof path, "/backends/%s/status", name);
        if (get_checked(s, path, status_json) != 0) return -1;
    }
    return 0;
}

/* ---- calibration parsing ---- */

static double num_of(cJSON *o, const char *key, double dflt)
{
    cJSON *v = cJSON_GetObjectItem(o, key);
    return cJSON_IsNumber(v) ? v->valuedouble : dflt;
}

/* Qubit index from a JSON value: real backends send numbers, the OpenAPI schema says strings. */
static int json_qubit(const cJSON *v)
{
    if (cJSON_IsNumber(v)) return v->valueint;
    if (cJSON_IsString(v) && v->valuestring) return atoi(v->valuestring);
    return -1;
}

static int find_edge(const IBMCalibration *cal, int a, int b)
{
    for (int e = 0; e < cal->nedge; e++)
        if ((cal->edge_a[e] == a && cal->edge_b[e] == b) || (cal->edge_a[e] == b && cal->edge_b[e] == a))
            return e;
    return -1;
}

int ibm_calibration_parse(const char *cfg_json, const char *props_json, IBMCalibration *cal)
{
    memset(cal, 0, sizeof *cal);
    for (int q = 0; q < IBM_MAX_QUBIT; q++) {
        cal->readout_err[q] = cal->sx_err[q] = cal->x_err[q] = NAN;
        cal->sx_len_ns[q] = NAN;
        cal->t1_us[q] = cal->t2_us[q] = NAN;
    }
    cJSON *cfg = cJSON_Parse(cfg_json);
    cJSON *pr = cJSON_Parse(props_json);
    if (!cfg || !pr) {
        cJSON_Delete(cfg);
        cJSON_Delete(pr);
        return -1;
    }
    cal->nqubit = (int)num_of(cfg, "n_qubits", 0);
    if (cal->nqubit <= 0 || cal->nqubit > IBM_MAX_QUBIT) {
        cJSON_Delete(cfg);
        cJSON_Delete(pr);
        return -1;
    }
    cJSON *bn = cJSON_GetObjectItem(cfg, "backend_name");
    cJSON *bv = cJSON_GetObjectItem(cfg, "backend_version");
    snprintf(cal->backend_name, sizeof cal->backend_name, "%s", cJSON_IsString(bn) ? bn->valuestring : "");
    snprintf(cal->backend_version, sizeof cal->backend_version, "%s", cJSON_IsString(bv) ? bv->valuestring : "");

    cJSON *basis = cJSON_GetObjectItem(cfg, "basis_gates");
    for (cJSON *g = basis ? basis->child : NULL; g; g = g->next) {
        if (cJSON_IsString(g) && strcmp(g->valuestring, "cz") == 0) cal->has_cz = 1;
        if (cJSON_IsString(g) && strcmp(g->valuestring, "ecr") == 0) cal->has_ecr = 1;
    }

    /* couplers: undirected unique pairs from the coupling map */
    cJSON *cm = cJSON_GetObjectItem(cfg, "coupling_map");
    for (cJSON *pair = cm ? cm->child : NULL; pair; pair = pair->next) {
        if (!cJSON_IsArray(pair) || cJSON_GetArraySize(pair) != 2) continue;
        int a = cJSON_GetArrayItem(pair, 0)->valueint;
        int b = cJSON_GetArrayItem(pair, 1)->valueint;
        if (a == b || a < 0 || b < 0 || a >= cal->nqubit || b >= cal->nqubit) continue;
        if (find_edge(cal, a, b) >= 0) continue;
        if (cal->nedge >= IBM_MAX_EDGE) break;
        cal->edge_a[cal->nedge] = a < b ? a : b;
        cal->edge_b[cal->nedge] = a < b ? b : a;
        cal->edge_err[cal->nedge] = NAN;
        cal->edge_len_ns[cal->nedge] = NAN;
        cal->nedge++;
    }

    /* qubit properties: qubits[q] = [ {name,value,unit,date}, ... ] */
    cJSON *qs = cJSON_GetObjectItem(pr, "qubits");
    int q = 0;
    for (cJSON *arr = qs ? qs->child : NULL; arr && q < cal->nqubit; arr = arr->next, q++) {
        for (cJSON *p = arr->child; p; p = p->next) {
            cJSON *nm = cJSON_GetObjectItem(p, "name");
            if (!cJSON_IsString(nm)) continue;
            double v = num_of(p, "value", NAN);
            if (strcmp(nm->valuestring, "readout_error") == 0) cal->readout_err[q] = v;
            else if (strcmp(nm->valuestring, "T1") == 0) cal->t1_us[q] = v;
            else if (strcmp(nm->valuestring, "T2") == 0) cal->t2_us[q] = v;
        }
    }

    /* gate properties: gates[] = {gate, qubits:[..], parameters:[{name,value,unit}]} */
    cJSON *gates = cJSON_GetObjectItem(pr, "gates");
    for (cJSON *g = gates ? gates->child : NULL; g; g = g->next) {
        cJSON *gn = cJSON_GetObjectItem(g, "gate");
        cJSON *gq = cJSON_GetObjectItem(g, "qubits");
        cJSON *gp = cJSON_GetObjectItem(g, "parameters");
        if (!cJSON_IsString(gn) || !cJSON_IsArray(gq)) continue;
        int nq = cJSON_GetArraySize(gq);
        int qa = nq > 0 ? json_qubit(cJSON_GetArrayItem(gq, 0)) : -1;
        int qb = nq > 1 ? json_qubit(cJSON_GetArrayItem(gq, 1)) : -1;
        double err = NAN, len = NAN;
        for (cJSON *p = gp ? gp->child : NULL; p; p = p->next) {
            cJSON *pn = cJSON_GetObjectItem(p, "name");
            if (!cJSON_IsString(pn)) continue;
            if (strcmp(pn->valuestring, "gate_error") == 0) err = num_of(p, "value", NAN);
            if (strcmp(pn->valuestring, "gate_length") == 0) len = num_of(p, "value", NAN);
        }
        if (strcmp(gn->valuestring, "cz") == 0 && qa >= 0 && qb >= 0) {
            int e = find_edge(cal, qa, qb);
            if (e >= 0) { cal->edge_err[e] = err; cal->edge_len_ns[e] = len; }
        } else if ((strcmp(gn->valuestring, "sx") == 0 || strcmp(gn->valuestring, "x") == 0) && qa >= 0 && qa < IBM_MAX_QUBIT) {
            if (strcmp(gn->valuestring, "sx") == 0) { cal->sx_err[qa] = err; cal->sx_len_ns[qa] = len; }
            else cal->x_err[qa] = err;
        }
    }
    cJSON_Delete(cfg);
    cJSON_Delete(pr);
    return 0;
}

/* ---- jobs ---- */

/* Estimator V2 job with several PUBs (one per circuit), all in one job so they share the
   backend, calibration window and options. labels[p] holds the observables of PUB p. */
int ibm_job_submit_pubs(IBMSession *s, const char *backend, const char *const *qasms, int npubs,
                        const char *const *const *labels, int nlabels, double precision,
                        int resilience_level, int seed, char **resp)
{
    cJSON *root = cJSON_CreateObject();
    cJSON_AddStringToObject(root, "program_id", "estimator");
    cJSON_AddStringToObject(root, "backend", backend);
    cJSON *params = cJSON_AddObjectToObject(root, "params");
    cJSON *pubs = cJSON_AddArrayToObject(params, "pubs");
    for (int p = 0; p < npubs; p++) {
        cJSON *pub = cJSON_CreateArray();
        cJSON_AddItemToArray(pub, cJSON_CreateString(qasms[p]));
        cJSON *obs = cJSON_CreateArray();
        for (int i = 0; i < nlabels; i++) cJSON_AddItemToArray(obs, cJSON_CreateString(labels[p][i]));
        cJSON_AddItemToArray(pub, obs);
        cJSON_AddItemToArray(pubs, pub);
    }
    cJSON_AddNumberToObject(params, "version", 2);
    cJSON_AddNumberToObject(params, "resilience_level", resilience_level);
    cJSON *opts = cJSON_AddObjectToObject(params, "options");
    cJSON_AddNumberToObject(opts, "default_precision", precision);
    if (seed >= 0) cJSON_AddNumberToObject(opts, "seed_estimator", seed);
    char *body = cJSON_PrintUnformatted(root);
    cJSON_Delete(root);
    long code = 0;
    char *r = NULL;
    int rc = ibm_http(s, "POST", "/jobs", body, &r, &code);
    free(body);
    if (rc != 0 || (code != 200 && code != 201)) {
        ibm_report_http_error(s, "POST", "/jobs", code, r);
        if (resp) *resp = r; else free(r);
        return -1;
    }
    if (resp) *resp = r; else free(r);
    return 0;
}

int ibm_job_submit(IBMSession *s, const char *backend, const char *qasm, const char *const *labels, int nlabels,
                   double precision, char **resp)
{
    const char *const *lab[1] = {labels};
    return ibm_job_submit_pubs(s, backend, &qasm, 1, lab, nlabels, precision, 1, -1, resp);
}

int ibm_job_get(IBMSession *s, const char *job_id, char **json)
{
    char path[256];
    long code = 0;
    snprintf(path, sizeof path, "/jobs/%s", job_id);
    (void)code;
    return get_checked(s, path, json);
}

int ibm_job_results(IBMSession *s, const char *job_id, char **json)
{
    char path[256];
    long code = 0;
    snprintf(path, sizeof path, "/jobs/%s/results", job_id);
    (void)code;
    return get_checked(s, path, json);
}

/* ---- base64 (little-endian f64 arrays) ---- */

static int b64val(int c)
{
    if (c >= 'A' && c <= 'Z') return c - 'A';
    if (c >= 'a' && c <= 'z') return c - 'a' + 26;
    if (c >= '0' && c <= '9') return c - '0' + 52;
    if (c == '+') return 62;
    if (c == '/') return 63;
    return -1;
}

int base64_decode_f64(const char *b64, double **out, size_t *n)
{
    size_t len = strlen(b64);
    unsigned char *bytes = qc_calloc(len, 1);
    size_t nb = 0;
    unsigned acc = 0;
    int bits = 0;
    for (size_t i = 0; i < len; i++) {
        if (b64[i] == '=') break;
        int v = b64val((unsigned char)b64[i]);
        if (v < 0) continue;
        acc = (acc << 6) | (unsigned)v;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            bytes[nb++] = (unsigned char)((acc >> bits) & 0xff);
        }
    }
    if (nb % 8 != 0) { free(bytes); return -1; }
    *n = nb / 8;
    *out = qc_calloc(*n ? *n : 1, sizeof(double));
    for (size_t k = 0; k < *n; k++) {
        uint64_t u = 0;
        for (int b = 7; b >= 0; b--) u = (u << 8) | bytes[k * 8 + (size_t)b];
        double d;
        memcpy(&d, &u, sizeof d);
        (*out)[k] = d;
    }
    free(bytes);
    return 0;
}
