/*
 * IBM credentials: entered via the CLI, stored in a user-only file.
 *
 *   qchem-vqe auth login                      prompts for API key (hidden) and CRN
 *   qchem-vqe auth login --api-key-stdin --crn CRN   key read from stdin (scripts)
 *   qchem-vqe auth status                     shows where credentials come from (key masked)
 *   qchem-vqe auth logout                     deletes the stored credentials file
 *
 * Storage: $QCHEM_VQE_CREDENTIALS, else $XDG_CONFIG_HOME/qchem-vqe/ibm_credentials,
 * else $HOME/.config/qchem-vqe/ibm_credentials. The directory is created 0700 and the
 * file 0600. Environment variables take precedence over the file.
 *
 * The API key is never accepted as a command-line argument (it would show up in
 * shell history and in `ps`) and is never printed in full.
 */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <termios.h>
#include <unistd.h>

#define CRED_FILE_NAME "ibm_credentials"

int cred_default_path(char *buf, size_t n)
{
    if (qc_env_get("QCHEM_VQE_CREDENTIALS", buf, n)) return 0;
    char base[768];
    if (qc_env_get("XDG_CONFIG_HOME", base, sizeof base)) {
        snprintf(buf, n, "%s/qchem-vqe/%s", base, CRED_FILE_NAME);
        return 0;
    }
    if (qc_env_get("HOME", base, sizeof base)) {
        snprintf(buf, n, "%s/.config/qchem-vqe/%s", base, CRED_FILE_NAME);
        return 0;
    }
    return -1;
}

static void trim(char *s)
{
    size_t n = strlen(s);
    while (n && (s[n - 1] == '\n' || s[n - 1] == '\r' || s[n - 1] == ' ' || s[n - 1] == '\t')) s[--n] = 0;
    char *p = s;
    while (*p == ' ' || *p == '\t') p++;
    if (p != s) memmove(s, p, strlen(p) + 1);
}

/* Reads the credentials file into s (only the keys present). Returns 0 if the file exists. */
int cred_load(const char *path, IBMSession *s)
{
    size_t len = 0;
    char *txt = qc_read_file(path, &len);
    if (!txt) return -1;
    char *save = NULL;
    for (char *line = strtok_r(txt, "\n", &save); line; line = strtok_r(NULL, "\n", &save)) {
        if (line[0] == '#' || line[0] == 0) continue;
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = 0;
        char *key = line, *val = eq + 1;
        trim(key);
        trim(val);
        if (strcmp(key, "IBM_QUANTUM_API_KEY") == 0) snprintf(s->api_key, sizeof s->api_key, "%s", val);
        else if (strcmp(key, "IBM_QUANTUM_CRN") == 0) snprintf(s->crn, sizeof s->crn, "%s", val);
        else if (strcmp(key, "IBM_QUANTUM_API_URL") == 0) snprintf(s->base_url, sizeof s->base_url, "%s", val);
        else if (strcmp(key, "IBM_QUANTUM_API_VERSION") == 0) snprintf(s->api_version, sizeof s->api_version, "%s", val);
    }
    free(txt);
    return 0;
}

/* Writes the credentials file with user-only permissions. */
int cred_store(const char *path, const char *api_key, const char *crn, const char *url, const char *version)
{
    char dir[1024];
    snprintf(dir, sizeof dir, "%s", path);
    char *slash = strrchr(dir, '/');
    if (slash) {
        *slash = 0;
        if (qc_mkdir_p(dir) != 0) return -1;
        chmod(dir, 0700);
    }
    int fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (fd < 0) return -1;
    fchmod(fd, 0600);
    FILE *f = fdopen(fd, "w");
    if (!f) { close(fd); return -1; }
    fprintf(f, "# qchem-vqe IBM Quantum credentials (written by `qchem-vqe auth login`)\n");
    fprintf(f, "IBM_QUANTUM_API_KEY=%s\n", api_key);
    fprintf(f, "IBM_QUANTUM_CRN=%s\n", crn);
    if (url && *url) fprintf(f, "IBM_QUANTUM_API_URL=%s\n", url);
    if (version && *version) fprintf(f, "IBM_QUANTUM_API_VERSION=%s\n", version);
    int ok = fclose(f) == 0;
    return ok ? 0 : -1;
}

/* "****abcd": the last four characters only */
void cred_mask(const char *secret, char *out, size_t n)
{
    size_t len = strlen(secret);
    if (len <= 4) snprintf(out, n, "****");
    else snprintf(out, n, "****%s", secret + len - 4);
}

/* Reads a line; if echo can be disabled (terminal), the input is not shown. */
static int read_line(const char *label, int hidden, char *buf, size_t n)
{
    fprintf(stderr, "%s", label);
    fflush(stderr);
    struct termios old;
    int tty = isatty(STDIN_FILENO);
    int changed = 0;
    if (hidden && tty && tcgetattr(STDIN_FILENO, &old) == 0) {
        struct termios t = old;
        t.c_lflag &= ~ECHO;
        if (tcsetattr(STDIN_FILENO, TCSANOW, &t) == 0) changed = 1;
    }
    char *r = fgets(buf, (int)n, stdin);
    if (changed) {
        tcsetattr(STDIN_FILENO, TCSANOW, &old);
        fprintf(stderr, "\n");
    }
    if (!r) return -1;
    trim(buf);
    return 0;
}

static void usage_auth(void)
{
    fprintf(stderr,
            "usage:\n"
            "  qchem-vqe auth login [--crn CRN] [--url URL] [--no-verify]\n"
            "  qchem-vqe auth login --api-key-stdin --crn CRN [--url URL] [--no-verify]\n"
            "  qchem-vqe auth status\n"
            "  qchem-vqe auth logout\n");
}

static int cmd_login(int argc, char **argv)
{
    char crn[512] = "", url[256] = "", key[256] = "";
    int from_stdin = 0, verify = 1;
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--crn") == 0 && i + 1 < argc) snprintf(crn, sizeof crn, "%s", argv[++i]);
        else if (strcmp(argv[i], "--url") == 0 && i + 1 < argc) snprintf(url, sizeof url, "%s", argv[++i]);
        else if (strcmp(argv[i], "--api-key-stdin") == 0) from_stdin = 1;
        else if (strcmp(argv[i], "--no-verify") == 0) verify = 0;
        else if (strcmp(argv[i], "--api-key") == 0) {
            fprintf(stderr, "error: --api-key is refused: a key on the command line leaks into shell history. "
                            "Run `auth login` (hidden prompt) or pipe it with --api-key-stdin.\n");
            return 1;
        } else { usage_auth(); return 1; }
    }

    if (from_stdin) {
        if (!fgets(key, sizeof key, stdin)) { fprintf(stderr, "error: no API key on stdin\n"); return 1; }
        trim(key);
        if (!crn[0]) { fprintf(stderr, "error: --api-key-stdin needs --crn CRN\n"); return 1; }
    } else {
        if (!isatty(STDIN_FILENO)) {
            fprintf(stderr, "error: stdin is not a terminal; use --api-key-stdin for non-interactive login\n");
            return 1;
        }
        if (read_line("IBM Cloud API key (input hidden): ", 1, key, sizeof key) != 0 || !key[0]) {
            fprintf(stderr, "error: no API key entered\n");
            return 1;
        }
        if (!crn[0]) {
            if (read_line("Quantum instance CRN: ", 0, crn, sizeof crn) != 0 || !crn[0]) {
                fprintf(stderr, "error: no CRN entered\n");
                return 1;
            }
        }
    }
    if (strlen(key) < 8) { fprintf(stderr, "error: API key looks too short\n"); return 1; }

    IBMSession s;
    memset(&s, 0, sizeof s);
    snprintf(s.api_key, sizeof s.api_key, "%s", key);
    snprintf(s.crn, sizeof s.crn, "%s", crn);
    if (verify) {
        printf("verifying the API key with IBM Cloud IAM...\n");
        if (ibm_session_authenticate(&s) != 0) {
            fprintf(stderr, "error: IAM rejected the key; nothing was saved (use --no-verify to save anyway)\n");
            memset(key, 0, sizeof key);
            memset(s.api_key, 0, sizeof s.api_key);
            return 1;
        }
        printf("  key accepted\n");
    }

    char path[1024];
    if (cred_default_path(path, sizeof path) != 0) { fprintf(stderr, "error: no HOME or XDG_CONFIG_HOME\n"); return 1; }
    if (cred_store(path, key, crn, url, NULL) != 0) { fprintf(stderr, "error: cannot write %s\n", path); return 1; }
    char masked[32];
    cred_mask(key, masked, sizeof masked);
    memset(key, 0, sizeof key);
    memset(s.api_key, 0, sizeof s.api_key);
    printf("saved credentials to %s (mode 0600)\n", path);
    printf("  API key: %s   CRN: %s\n", masked, crn);
    return 0;
}

static int cmd_status(void)
{
    char path[1024];
    cred_default_path(path, sizeof path);
    IBMSession file;
    memset(&file, 0, sizeof file);
    int have_file = cred_load(path, &file) == 0;
    char k[256] = "", c[512] = "";
    int env_key = qc_env_get("IBM_QUANTUM_API_KEY", k, sizeof k);
    int env_crn = qc_env_get("IBM_QUANTUM_CRN", c, sizeof c);
    char masked[32];
    printf("credentials file : %s (%s)\n", path, have_file ? "present" : "absent");
    if (have_file) {
        cred_mask(file.api_key, masked, sizeof masked);
        printf("  API key        : %s\n", file.api_key[0] ? masked : "(not set)");
        printf("  CRN            : %s\n", file.crn[0] ? file.crn : "(not set)");
    }
    printf("environment      : IBM_QUANTUM_API_KEY %s, IBM_QUANTUM_CRN %s\n",
           env_key ? "set" : "unset", env_crn ? "set" : "unset");
    if (env_key) {
        cred_mask(k, masked, sizeof masked);
        printf("  env API key    : %s (overrides the file)\n", masked);
    }
    memset(k, 0, sizeof k);
    int usable = (file.api_key[0] && file.crn[0]) || (env_key && env_crn);
    if (usable) {
        IBMSession eff;
        if (ibm_session_env(&eff) == 0) {
            char region[64] = "?";
            ibm_crn_region(eff.crn, region, sizeof region);
            printf("instance region  : %s\n", region);
            printf("API endpoint     : %s\n", eff.base_url);
            memset(eff.api_key, 0, sizeof eff.api_key);
        }
    }
    printf("status           : %s\n", usable ? "credentials available" : "no credentials (run `qchem-vqe auth login`)");
    return 0;
}

static int cmd_logout(void)
{
    char path[1024];
    if (cred_default_path(path, sizeof path) != 0) return 1;
    if (remove(path) == 0) printf("removed %s\n", path);
    else if (errno == ENOENT) printf("no credentials file at %s\n", path);
    else { fprintf(stderr, "error: cannot remove %s\n", path); return 1; }
    return 0;
}

int cmd_auth(int argc, char **argv)
{
    if (argc < 1) { usage_auth(); return 1; }
    if (strcmp(argv[0], "login") == 0) return cmd_login(argc - 1, argv + 1);
    if (strcmp(argv[0], "status") == 0) return cmd_status();
    if (strcmp(argv[0], "logout") == 0) return cmd_logout();
    usage_auth();
    return 1;
}
