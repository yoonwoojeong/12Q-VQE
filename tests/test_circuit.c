/*
 * Simulator, circuit IR, native decomposition, UCCSD ansatz checks.
 */
#include "testlib.h"

static double overlap_mag(const StateVec *a, const StateVec *b)
{
    double complex s = 0;
    for (size_t i = 0; i < a->dim; i++) s += conj(a->a[i]) * b->a[i];
    return cabs(s);
}

int main(void)
{
    printf("== test_circuit\n");

    /* 1. Pauli rotations: unitarity and comparison with gate-level compilation */
    StateVec a, b, c;
    sv_init(&a, 4); sv_init(&b, 4); sv_init(&c, 4);
    /* build a non-trivial state */
    double complex h[2][2] = {{0.70710678118654752, 0.70710678118654752}, {0.70710678118654752, -0.70710678118654752}};
    for (int q = 0; q < 4; q++) sv_apply_1q(&a, q, h);
    sv_apply_cx(&a, 0, 2);
    sv_apply_cz(&a, 1, 3);
    sv_copy(&b, &a); sv_copy(&c, &a);
    uint32_t xs[] = {0x3, 0x5, 0xf, 0x9};
    uint32_t zs[] = {0x0, 0xa, 0x6, 0x1};
    double phis[] = {0.31, -1.2, 0.77, 2.4};
    for (int k = 0; k < 4; k++) {
        sv_apply_pauli_rotation(&a, xs[k], zs[k], phis[k]);
    }
    CHECK(fabs(sv_norm(&a) - 1.0) < 1e-12, "Pauli rotations preserve the norm (%.15f)", sv_norm(&a));

    Circuit ic;
    circ_init(&ic, 4);
    for (int k = 0; k < 4; k++) circ_append_pauli_rotation(&ic, xs[k], zs[k], phis[k], 4);
    circ_simulate(&ic, &b);
    /* the IR was built from the same state prep, so compare against c evolved via IR */
    Circuit nat;
    circ_to_native(&ic, &nat);
    circ_simulate(&nat, &c);
    /* apply IR to the prepared state for each case: rebuild state prep on both sides */
    StateVec p1, p2;
    sv_init(&p1, 4); sv_init(&p2, 4);
    for (int q = 0; q < 4; q++) sv_apply_1q(&p1, q, h);
    sv_apply_cx(&p1, 0, 2); sv_apply_cz(&p1, 1, 3);
    sv_copy(&p2, &p1);
    circ_simulate(&ic, &p1);          /* gate IR of the four rotations */
    circ_simulate(&nat, &p2);         /* native-gate decomposition of the same */
    /* direct Pauli-exponential evolution of the same prep */
    StateVec d;
    sv_init(&d, 4);
    for (int q = 0; q < 4; q++) sv_apply_1q(&d, q, h);
    sv_apply_cx(&d, 0, 2); sv_apply_cz(&d, 1, 3);
    for (int k = 0; k < 4; k++) sv_apply_pauli_rotation(&d, xs[k], zs[k], phis[k]);
    double ov_ir = overlap_mag(&d, &p1);
    double ov_nat = overlap_mag(&d, &p2);
    CHECK(fabs(ov_ir - 1.0) < 1e-10, "gate IR of Pauli rotations == direct exp(iphiP) (|<d|ir>| = %.12f)", ov_ir);
    CHECK(fabs(ov_nat - 1.0) < 1e-10, "native rz/sx/x/cz decomposition == direct (up to global phase) (%.12f)", ov_nat);
    CHECK(circ_count_kind(&nat, G_CX) == 0 && circ_count_kind(&nat, G_H) == 0,
          "native circuit contains only rz, sx, x, cz (CX/H removed)");
    (void)b; (void)c;

    /* 2. QASM 3 emission */
    char *q3 = qasm3_string(&nat, 4, 1);
    CHECK(strncmp(q3, "OPENQASM 3.0;", 13) == 0 && strstr(q3, "qubit[4] q;") != NULL && strstr(q3, "cz q[") != NULL,
          "OpenQASM 3 header and native gates emitted");
    free(q3);

    /* 3. UCCSD on CAS(2,6) */
    Sys s;
    sys_build(&s, 2, 6);
    UCCSD u;
    uccsd_build(&u, s.fh.nspin, 2);
    CHECK(u.nexc == 35, "UCCSD has 35 spin-conserving singles+doubles (got %d)", u.nexc);
    int nsingle = 0, ndouble = 0;
    for (int k = 0; k < u.nexc; k++) (u.exc[k].kind == 1 ? nsingle++ : ndouble++);
    CHECK(nsingle == 10 && ndouble == 25, "singles %d, doubles %d", nsingle, ndouble);
    double *th = qc_calloc((size_t)u.nexc, sizeof(double));
    StateVec sc;
    sv_init(&sc, s.fh.nspin);
    double e0 = vqe_energy(&u, &s.fh.ham, th, &sc);
    CHECK(fabs(e0 - s.rhf.energy) < 1e-10, "UCCSD at theta=0 is the HF energy (%.12f)", e0);

    FCIResult f;
    fci_solve(&f, &s.fh, &s.as, 2);
    unsigned seed = 3;
    int viol = 0;
    double emin = 1e9;
    for (int t = 0; t < 20; t++) {
        for (int k = 0; k < u.nexc; k++) {
            seed = seed * 1103515245u + 12345u;
            th[k] = (((seed >> 8) & 0xffff) / 65535.0 - 0.5) * 1.5;
        }
        double e = vqe_energy(&u, &s.fh.ham, th, &sc);
        if (e < f.ground - 1e-10) viol++;
        if (e < emin) emin = e;
    }
    CHECK(viol == 0, "variational principle: 20 random UCCSD states all >= FCI (min %.10f, FCI %.10f)", emin, f.ground);

    /* gate-level UCCSD circuit reproduces the direct state preparation */
    for (int k = 0; k < u.nexc; k++) th[k] = 0.05 * (k % 7) - 0.1;
    Circuit uc;
    uccsd_circuit(&u, th, &uc);
    StateVec s1, s2;
    sv_init(&s1, s.fh.nspin); sv_init(&s2, s.fh.nspin);
    circ_simulate(&uc, &s1);
    uccsd_prepare(&u, th, &s2);
    double ov = overlap_mag(&s1, &s2);
    CHECK(fabs(ov - 1.0) < 1e-10, "UCCSD gate circuit == direct Pauli-rotation preparation (|<.|.>| = %.12f)", ov);
    Circuit un;
    circ_to_native(&uc, &un);
    StateVec s3;
    sv_init(&s3, s.fh.nspin);
    circ_simulate(&un, &s3);
    ov = overlap_mag(&s1, &s3);
    CHECK(fabs(ov - 1.0) < 1e-9, "UCCSD native circuit == IR (|<.|.>| = %.12f); gates %d -> %d, 2q %d -> %d",
          ov, uc.ng, un.ng, circ_count_2q(&uc), circ_count_2q(&un));

    /* 4. Givens state preparation reproduces a two-electron state exactly */
    StateVec tg;
    sv_init(&tg, s.fh.nspin);
    memset(tg.a, 0, sizeof(double complex) * tg.dim);
    for (int i = 0; i < f.ndet; i++) tg.a[f.dets[i]] = f.ground_vec[i];
    Circuit gc;
    GivensInfo gi;
    givens_state_circuit(&tg, s.fh.nspin, 0.0, &gc, &gi);
    StateVec gv;
    sv_init(&gv, s.fh.nspin);
    circ_simulate(&gc, &gv);
    double eg = sv_expect_sum(&gv, &s.fh.ham);
    CHECK(gi.fidelity > 1.0 - 1e-10, "Givens circuit prepares the FCI state (fidelity %.12f)", gi.fidelity);
    CHECK(fabs(eg - f.ground) < 1e-9, "Givens circuit energy equals FCI (%.12f)", eg);
    CHECK(circ_count_kind(&gc, G_CX) < circ_count_2q(&uc) / 3,
          "Givens circuit is shallow: %d CX vs %d for UCCSD", circ_count_kind(&gc, G_CX), circ_count_2q(&uc));
    Circuit gt;
    GivensInfo gi2;
    givens_state_circuit(&tg, s.fh.nspin, 1e-4, &gt, &gi2);
    CHECK(gi2.fidelity > 0.9999 && gi2.r < 6, "truncated Givens circuit (r=%d) keeps fidelity %.8f", gi2.r, gi2.fidelity);
    circ_free(&gc);
    circ_free(&gt);
    sv_free(&tg);
    sv_free(&gv);

    /* 5. Analytic adjoint gradient matches central finite differences & ADAPT-VQE reaches FCI */
    {
        double *grad = qc_calloc((size_t)u.nexc, sizeof(double));
        StateVec psi_g, lam_g;
        sv_init(&psi_g, s.fh.nspin);
        sv_init(&lam_g, s.fh.nspin);
        vqe_energy_grad(&u, &s.fh.ham, th, grad, &psi_g, &lam_g);
        double max_gerr = 0.0;
        double eps = 1e-6;
        for (int k = 0; k < u.nexc; k++) {
            double orig = th[k];
            th[k] = orig + eps;
            double ep = vqe_energy(&u, &s.fh.ham, th, &sc);
            th[k] = orig - eps;
            double em = vqe_energy(&u, &s.fh.ham, th, &sc);
            th[k] = orig;
            double fd = (ep - em) / (2.0 * eps);
            double err = fabs(grad[k] - fd);
            if (err > max_gerr) max_gerr = err;
        }
        CHECK(max_gerr < 1e-8, "reverse-mode adjoint gradient matches finite difference (max |diff| %.2e)", max_gerr);
        free(grad);
        sv_free(&psi_g);
        sv_free(&lam_g);

        AdaptResult ar;
        Circuit ac;
        adapt_vqe_optimize(&u, &s.fh.ham, 15, 1e-6, 0, &ar, &ac);
        CHECK(fabs(ar.energy - f.ground) < 1e-10 && ar.n_ops <= 10,
              "ADAPT-VQE + L-BFGS reaches FCI in %d ops (|E-FCI| = %.2e, %d CX)",
              ar.n_ops, fabs(ar.energy - f.ground), circ_count_2q(&ac));
        adapt_result_free(&ar);
        circ_free(&ac);
    }

    fci_free(&f);
    uccsd_free(&u);
    circ_free(&uc); circ_free(&un); circ_free(&ic); circ_free(&nat);
    sv_free(&a); sv_free(&b); sv_free(&c); sv_free(&p1); sv_free(&p2); sv_free(&d);
    sv_free(&sc); sv_free(&s1); sv_free(&s2); sv_free(&s3);
    free(th);
    fermion_free(&s.fh);
    active_space_free(&s.as);
    rhf_free(&s.rhf);
    ao_integrals_free(&s.ints);
    basis_free(&s.bs);
    printf("%s\n", g_fail ? "test_circuit: FAILED" : "test_circuit: passed");
    return g_fail ? 1 : 0;
}
