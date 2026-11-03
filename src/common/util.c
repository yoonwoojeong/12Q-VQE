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
