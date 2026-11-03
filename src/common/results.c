/* Result files: JSON with experiment metadata. */
#define _POSIX_C_SOURCE 200809L
#include "qchem.h"

#include <cjson/cJSON.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef QC_GIT_REV
#define QC_GIT_REV "unknown"
#endif

int results_write_json(const char *path, const char *json_text)
{
    char dir[1024];
    snprintf(dir, sizeof dir, "%s", path);
    char *slash = strrchr(dir, '/');
    if (slash) {
        *slash = 0;
        qc_mkdir_p(dir);
    }
    size_t n = strlen(json_text);
    if (qc_write_file(path, json_text, n) != 0) {
        qc_warn("could not write '%s'", path);
        return -1;
    }
    printf("  wrote %s\n", path);
    return 0;
}

char *results_metadata_json(void)
{
    cJSON *m = cJSON_CreateObject();
    cJSON_AddStringToObject(m, "program", "qchem-vqe");
    cJSON_AddStringToObject(m, "version", QC_VERSION);
    cJSON_AddStringToObject(m, "git_revision", QC_GIT_REV);
    char ts[64];
    time_t now = time(NULL);
    struct tm tmv;
    gmtime_r(&now, &tmv);
    strftime(ts, sizeof ts, "%Y-%m-%dT%H:%M:%SZ", &tmv);
    cJSON_AddStringToObject(m, "timestamp_utc", ts);
    char *s = cJSON_PrintUnformatted(m);
    cJSON_Delete(m);
    /* return a copy the caller can free with free() */
    char *out = qc_strdup(s);
    cJSON_free(s);
    return out;
}
