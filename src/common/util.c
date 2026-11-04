/* Common helpers: errors, memory, files, time, environment. */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>

void qc_die(const char *fmt, ...)
{
    va_list ap;
    fprintf(stderr, "error: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    exit(1);
}

void qc_warn(const char *fmt, ...)
{
    va_list ap;
    fprintf(stderr, "warning: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
}

void *qc_calloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1, sz ? sz : 1);
    if (!p) qc_die("out of memory (%zu x %zu bytes)", n, sz);
    return p;
}

void *qc_realloc(void *p, size_t n)
{
    void *q = realloc(p, n ? n : 1);
    if (!q) qc_die("out of memory (realloc %zu bytes)", n);
    return q;
}

char *qc_strdup(const char *s)
{
    size_t n = strlen(s) + 1;
    char *d = qc_calloc(n, 1);
    memcpy(d, s, n);
    return d;
}

int qc_mkdir_p(const char *path)
{
    char tmp[1024];
    size_t len = strlen(path);
    if (len == 0 || len >= sizeof tmp) return -1;
    memcpy(tmp, path, len + 1);
    for (size_t i = 1; i < len; i++) {
        if (tmp[i] == '/') {
            tmp[i] = 0;
            if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
            tmp[i] = '/';
        }
    }
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
    return 0;
}

char *qc_read_file(const char *path, size_t *len_out)
{
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    rewind(f);
    char *buf = qc_calloc((size_t)n + 1, 1);
    size_t got = fread(buf, 1, (size_t)n, f);
    fclose(f);
    buf[got] = 0;
    if (len_out) *len_out = got;
    return buf;
}

int qc_write_file(const char *path, const char *data, size_t len)
{
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    size_t w = fwrite(data, 1, len, f);
    fclose(f);
    return w == len ? 0 : -1;
}

double qc_wall_time(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int qc_env_get(const char *name, char *buf, size_t buflen)
{
    const char *v = getenv(name);
    if (!v || !*v) return 0;
    snprintf(buf, buflen, "%s", v);
    return 1;
}

void qc_print_banner(const char *what)
{
    printf("== qchem-vqe %s : %s ==\n", QC_VERSION, what);
}

static int dir_exists(const char *p)
{
    struct stat st;
    return stat(p, &st) == 0 && S_ISDIR(st.st_mode);
}

const char *qc_root_dir(char *buf, size_t n)
{
    if (qc_env_get("QCHEM_ROOT_DIR", buf, n)) return buf;
    char exe[1024] = {0};
    ssize_t len = readlink("/proc/self/exe", exe, sizeof(exe) - 1);
    if (len > 0) {
        exe[len] = 0;
        /* Check standard FHS install prefix: <prefix>/bin/qchem-vqe -> <prefix>/share/qchem-vqe/data/basis */
        char *bin_slash = strrchr(exe, '/');
        if (bin_slash) {
            *bin_slash = 0;
            char *prefix_slash = strrchr(exe, '/');
            if (prefix_slash && strcmp(prefix_slash, "/bin") == 0) {
                char share_probe[1200];
                snprintf(share_probe, sizeof share_probe, "%.*s/share/qchem-vqe/data/basis",
                         (int)(prefix_slash - exe), exe);
                if (dir_exists(share_probe)) {
                    snprintf(buf, n, "%.*s/share/qchem-vqe", (int)(prefix_slash - exe), exe);
                    return buf;
                }
            }
            *bin_slash = '/';
        }
        /* Walk up from executable path looking for data/basis */
        for (int up = 0; up < 4; up++) {
            char *slash = strrchr(exe, '/');
            if (!slash) break;
            *slash = 0;
            char probe[1200];
            snprintf(probe, sizeof probe, "%s/data/basis", exe);
            if (dir_exists(probe)) {
                snprintf(buf, n, "%s", exe);
                return buf;
            }
        }
    }
#ifdef QC_DATA_DIR
    {
        char tmp[1024];
        snprintf(tmp, sizeof tmp, "%s", QC_DATA_DIR);
        char *slash = strrchr(tmp, '/');
        if (slash && strcmp(slash, "/data") == 0) {
            *slash = 0;
            if (dir_exists(tmp)) {
                snprintf(buf, n, "%s", tmp);
                return buf;
            }
        }
    }
#endif
    snprintf(buf, n, ".");
    return buf;
}

const char *qc_data_dir(char *buf, size_t n)
{
    if (qc_env_get("QCHEM_DATA_DIR", buf, n)) return buf;
    char root[900];
    qc_root_dir(root, sizeof root);
    snprintf(buf, n, "%s/data", root);
    return buf;
}

const char *qc_cache_dir(char *buf, size_t n)
{
    if (qc_env_get("QCHEM_CACHE_DIR", buf, n)) return buf;
    char root[900];
    qc_root_dir(root, sizeof root);
    snprintf(buf, n, "%s/build/cache", root);
    return buf;
}

const char *qc_results_dir(char *buf, size_t n)
{
    if (qc_env_get("QCHEM_RESULTS_DIR", buf, n)) return buf;
    char root[900];
    qc_root_dir(root, sizeof root);
    snprintf(buf, n, "%s/results", root);
    return buf;
}
