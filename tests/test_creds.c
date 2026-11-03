#define _POSIX_C_SOURCE 200809L
/*
 * Credential storage: round trip, file permissions, masking, env precedence rules
 * (no network).
 */
#include "testlib.h"
#include <sys/stat.h>
#include <unistd.h>

int main(void)
{
    printf("== test_creds\n");
    char dir[] = "/tmp/qchem_creds_XXXXXX";
    CHECK(mkdtemp(dir) != NULL, "temporary directory created");
    char path[1024];
    snprintf(path, sizeof path, "%s/sub/ibm_credentials", dir);

    CHECK(cred_store(path, "KEY-ABCDEFGH-1234", "crn:v1:test::", "https://example.invalid/api/v1", "2026-04-15") == 0,
          "credentials stored");
    struct stat st;
    CHECK(stat(path, &st) == 0 && (st.st_mode & 0777) == 0600, "credentials file mode 0600");
    char dpath[1024];
    snprintf(dpath, sizeof dpath, "%s/sub", dir);
    CHECK(stat(dpath, &st) == 0 && (st.st_mode & 0777) == 0700, "credentials directory mode 0700");

    IBMSession s;
    memset(&s, 0, sizeof s);
    CHECK(cred_load(path, &s) == 0, "credentials loaded");
    CHECK(strcmp(s.api_key, "KEY-ABCDEFGH-1234") == 0, "API key round-trips");
    CHECK(strcmp(s.crn, "crn:v1:test::") == 0, "CRN round-trips");
    CHECK(strcmp(s.base_url, "https://example.invalid/api/v1") == 0, "API URL round-trips");
    CHECK(strcmp(s.api_version, "2026-04-15") == 0, "API version round-trips");

    char m[32];
    cred_mask("KEY-ABCDEFGH-1234", m, sizeof m);
    CHECK(strcmp(m, "****1234") == 0, "key masked to last four characters (%s)", m);
    cred_mask("abc", m, sizeof m);
    CHECK(strcmp(m, "****") == 0, "short secrets fully masked");

    char region[64];
    CHECK(ibm_crn_region("crn:v1:bluemix:public:quantum-computing:us-east:a/abc:def::", region, sizeof region) == 0 &&
          strcmp(region, "us-east") == 0, "CRN region parsed (us-east)");
    CHECK(strcmp(ibm_url_for_crn("crn:v1:bluemix:public:quantum-computing:us-east:a/abc:def::"),
                 "https://quantum.cloud.ibm.com/api/v1") == 0, "us-east instance uses the default endpoint");
    CHECK(strcmp(ibm_url_for_crn("crn:v1:bluemix:public:quantum-computing:eu-de:a/abc:def::"),
                 "https://eu-de.quantum.cloud.ibm.com/api/v1") == 0, "eu-de instance uses the eu-de endpoint");
    CHECK(ibm_crn_region("not-a-crn", region, sizeof region) != 0, "malformed CRN rejected");

    IBMSession none;
    memset(&none, 0, sizeof none);
    CHECK(cred_load("/nonexistent/qchem/ibm_credentials", &none) != 0, "missing file reported");

    unlink(path);
    rmdir(dpath);
    rmdir(dir);
    printf("%s\n", g_fail ? "test_creds: FAILED" : "test_creds: passed");
    return g_fail ? 1 : 0;
}
