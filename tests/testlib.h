/* Shared helpers for the scientific test programs. */
#ifndef TESTLIB_H
#define TESTLIB_H
#include "qchem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0;
#define CHECK(cond, ...) do { if (cond) { printf("  ok   : " __VA_ARGS__); printf("\n"); } \
    else { printf("  FAIL : " __VA_ARGS__); printf("\n"); g_fail++; } } while (0)

/* H2 at 0.735 A, cc-pVQZ, RHF, CAS(2,6) active space: builds everything */
typedef struct {
    Molecule mol;
    BasisSet bs;
    AOInts ints;
    RHFResult rhf;
    ActiveSpace as;
    FermionHam fh;
} Sys;

static inline void sys_build(Sys *s, int nel_act, int norb)
{
    memset(s, 0, sizeof *s);
    mol_h2_setup(&s->mol, 0.735, 0, 1);
    basis_build(&s->bs, &s->mol, QC_DATA_DIR "/basis/cc-pVQZ.nw", "cc-pVQZ");
    ao_integrals_compute(&s->ints, &s->bs, &s->mol);
    rhf_run(&s->rhf, &s->ints, &s->mol, 200, 1e-11, 1e-9, 0);
    int first = (s->mol.nelec - nel_act) / 2;
    active_space_build(&s->as, &s->ints, &s->rhf, &s->mol, nel_act, norb, first);
    fermion_build_hamiltonian(&s->fh, &s->as);
}

#endif
