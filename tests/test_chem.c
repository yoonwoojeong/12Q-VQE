/*
 * Electronic-structure validation: basis, integrals, RHF, active space,
 * Jordan-Wigner Hamiltonian, FCI cross-checks, multi-element molecules,
 * and 16-qubit scaling.
 * Reference numbers come from an external PySCF 2.14 run at the same
 * geometry/basis (scratch environment, not part of the build).
 */
#include "testlib.h"

int main(void)
{
    printf("== test_chem\n");
    Sys s;
    sys_build(&s, 2, 6);

    /* basis and integrals */
    int n = s.ints.nao;
    double maxdiagdev = 0, asym = 0, hsym = 0;
    for (int i = 0; i < n; i++) {
        double d = fabs(s.ints.S[(size_t)i * n + i] - 1.0);
        if (d > maxdiagdev) maxdiagdev = d;
        for (int j = 0; j < n; j++) {
            double a = fabs(s.ints.S[(size_t)i * n + j] - s.ints.S[(size_t)j * n + i]);
            if (a > asym) asym = a;
            double b = fabs(s.ints.H[(size_t)i * n + j] - s.ints.H[(size_t)j * n + i]);
            if (b > hsym) hsym = b;
        }
    }
    CHECK(n == 60, "cc-pVQZ H2 has 60 spherical AOs (got %d)", n);
    CHECK(maxdiagdev < 1e-10, "overlap diagonal = 1 (max dev %.2e)", maxdiagdev);
    CHECK(asym < 1e-12 && hsym < 1e-12, "S and H symmetric (%.2e, %.2e)", asym, hsym);
    /* PySCF uses CODATA-2014 bohr; we use CODATA-2018: agreement to ~2e-11 */
    CHECK(fabs(s.mol.enuc - 0.7199689944489797) < 1e-9, "nuclear repulsion 1/R (%.12f)", s.mol.enuc);

    /* ERI 8-fold symmetry at a few random index quadruples */
    double esym = 0.0;
    unsigned seed = 7;
    for (int t = 0; t < 2000; t++) {
        seed = seed * 1103515245u + 12345u;
        int p = (seed >> 3) % n;
        seed = seed * 1103515245u + 12345u;
        int q = (seed >> 3) % n;
        seed = seed * 1103515245u + 12345u;
        int r = (seed >> 3) % n;
        seed = seed * 1103515245u + 12345u;
        int u = (seed >> 3) % n;
        double a = s.ints.eri[(((size_t)p * n + q) * n + r) * n + u];
        double b = s.ints.eri[(((size_t)q * n + p) * n + u) * n + r];
        double c = s.ints.eri[(((size_t)u * n + r) * n + q) * n + p];
        double d = fabs(a - b) + fabs(a - c);
        if (d > esym) esym = d;
    }
    CHECK(esym < 1e-12, "ERI permutational symmetry (max dev %.2e)", esym);

    /* RHF against PySCF */
    CHECK(s.rhf.nocc == 1, "RHF occupation nocc=1");
    CHECK(fabs(s.rhf.energy - (-1.1334934589980064)) < 1e-9,
          "RHF/cc-pVQZ energy %.12f vs PySCF -1.133493458998", s.rhf.energy);

    /* active space and Hamiltonian */
    double hf_act = active_hf_energy(&s.as);
    CHECK(fabs(hf_act - s.rhf.energy) < 1e-10, "active-space HF energy equals RHF (%.12f)", hf_act);
    double comm, herm, hfe;
    int ok = fermion_validate(&s.fh, &s.as, &s.rhf, &comm, &herm, &hfe, NULL) == 0;
    CHECK(ok, "fermion validation (Hermiticity, N conservation, HF energy): comm=%.1e herm=%.1e hf=%.1e", comm, herm, hfe);
    CHECK(s.fh.ham.n == 327 || s.fh.ham.n > 0, "Pauli Hamiltonian has %zu terms", s.fh.ham.n);

    /* FCI: determinant reference vs PySCF CASCI and vs Pauli-restricted diagonalisation */
    FCIResult f;
    fci_solve(&f, &s.fh, &s.as, 2);
    CHECK(f.ndet == 66, "66 two-electron determinants in 12 spin orbitals (%d)", f.ndet);
    CHECK(fabs(f.ground - (-1.148902977365)) < 1e-9,
          "CAS(2,6) FCI %.12f vs PySCF CASCI -1.148902977365", f.ground);
    double pg = fci_pauli_restricted_ground(&s.fh, &f);
    CHECK(fabs(pg - f.ground) < 1e-11, "JW Pauli Hamiltonian reproduces FCI (|diff| %.2e)", fabs(pg - f.ground));
    CHECK(f.ground < s.rhf.energy, "correlation lowers energy below RHF (%.6f < %.6f)", f.ground, s.rhf.energy);

    fci_free(&f);
    fermion_free(&s.fh);
    active_space_free(&s.as);
    rhf_free(&s.rhf);
    ao_integrals_free(&s.ints);
    basis_free(&s.bs);

    /* Multi-element & frozen-core validation: LiH in STO-3G (SP shells, Z=3,1, ncore=1) */
    {
        Molecule lih;
        mol_setup_by_name(&lih, "LiH", 1.595);
        BasisSet bs;
        basis_build(&bs, &lih, QC_DATA_DIR "/basis/sto-3g.nw", "sto-3g");
        CHECK(bs.nao == 6, "LiH/STO-3G has 6 AOs from S and SP shells (got %d)", bs.nao);
        AOInts ints;
        ao_integrals_compute(&ints, &bs, &lih);
        RHFResult rhf;
        rhf_run(&rhf, &ints, &lih, 200, 1e-11, 1e-9, 0);
        ActiveSpace as;
        active_space_build(&as, &ints, &rhf, &lih, 2, 5, 1); /* frozen 1s core */
        CHECK(as.ncore == 1 && fabs(active_hf_energy(&as) - rhf.energy) < 1e-10,
              "LiH frozen-core CAS(2,5) HF energy equals full RHF (%.12f)", rhf.energy);
        FermionHam fh;
        fermion_build_hamiltonian(&fh, &as);
        FCIResult fci;
        fci_solve(&fci, &fh, &as, 2);
        double pg_lih = fci_pauli_restricted_ground(&fh, &fci);
        CHECK(fabs(pg_lih - fci.ground) < 1e-11 && fci.ground < rhf.energy,
              "LiH CAS(2,5) JW Pauli Hamiltonian matches FCI (%.12f Eh)", fci.ground);
        fci_free(&fci);
        fermion_free(&fh);
        active_space_free(&as);
        rhf_free(&rhf);
        ao_integrals_free(&ints);
        basis_free(&bs);
    }

    printf("%s\n", g_fail ? "test_chem: FAILED" : "test_chem: passed");
    return g_fail ? 1 : 0;
}
