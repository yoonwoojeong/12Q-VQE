# qchem-vqe

C17 implementation of a variational quantum eigensolver for H₂ (R = 0.735 Å) in the
cc-pVQZ basis, using a CAS(2,6) active space (12 qubits), a libcint integral engine,
an in-house 12-qubit statevector simulator, a UCCSD ansatz optimised with COBYLA
(NLopt), and a direct IBM Quantum Compute Service REST client that submits
Estimator V2 jobs with OpenQASM 3 circuits.

Nothing in this repository has been run on IBM hardware yet. See
[Status](#status-and-limitations).

## Build

Requirements: GCC (C17), `liblapacke-dev`, `libopenblas-dev`, `libnlopt-dev`,
`libcjson-dev`, `libcurl4-openssl-dev`, `zstd`, `unzip`, `curl`.

```bash
make fetch-libcint   # libcint 6.1.3 from conda-forge into third_party/libcint (not packaged by Ubuntu)
make                 # builds build/qchem-vqe
make test            # scientific validation tests
```

## Commands

```bash
./build/qchem-vqe prepare --basis cc-pVQZ          # libcint AO integrals, normalisation checks
./build/qchem-vqe rhf                              # RHF in the full cc-pVQZ basis (60 AO)
./build/qchem-vqe active-space --electrons 2 --orbitals 6 [--first N]
./build/qchem-vqe vqe --mode local [--maxeval N]   # UCCSD + COBYLA on the statevector simulator
./build/qchem-vqe hardware list
./build/qchem-vqe hardware inspect                 # backend selection, calibration summary, layout
./build/qchem-vqe hardware preflight               # routed circuit cost on the chosen layout
./build/qchem-vqe vqe --mode hardware              # dry run: builds and prints the plan
./build/qchem-vqe vqe --mode hardware --submit     # asks you to type SUBMIT; --yes skips the prompt
./build/qchem-vqe job status JOB_ID
./build/qchem-vqe job results JOB_ID
```

Results are written to `results/` (JSON with metadata: version, git revision, UTC time).
Hardware submissions are recorded in `results/hardware/`, including every observable
and its coefficient, so the job can be evaluated later without recomputation.

### Credentials

Read from the environment only, never printed or stored:

- `IBM_QUANTUM_API_KEY` — IBM Cloud API key (exchanged for an IAM bearer token)
- `IBM_QUANTUM_CRN` — Quantum instance CRN (`Service-CRN` header)
- `IBM_QUANTUM_API_URL` — optional, default `https://quantum.cloud.ibm.com/api/v1`
  (use `https://eu-de.quantum.cloud.ibm.com/api/v1` for the EU region)
- `IBM_QUANTUM_API_VERSION` — optional, default `2026-04-15`

## Results so far (local, no hardware)

| Quantity | Value (Eh) | Reference |
|---|---|---|
| Nuclear repulsion | 0.719968994426 | 1/R |
| RHF, full cc-pVQZ (60 AO) | −1.133493458998 | PySCF 2.14 RHF: −1.133493458998 |
| CAS(2,6) FCI (determinant, 66 dets) | −1.148902977328 | PySCF CASCI: −1.148902977365 |
| CAS(2,6) VQE, UCCSD (35 params), local | −1.148902977328 | equals FCI to 4×10⁻¹⁴ |
| CAS(2,3) FCI | −1.135444416491 | PySCF CASCI: −1.135444416507 |

The CAS(2,6) energy is an active-space quantity; the full cc-pVQZ treatment is the
RHF row. Both are reported separately in every output file.

## Method

- **Integrals**: libcint 6.1.3 spherical one- and two-electron integrals; cc-pVQZ for H
  from Basis Set Exchange (`data/basis/cc-pVQZ.nw`, `[4s3p2d1f]`, 30 functions per H).
  Each contracted shell is renormalised against libcint's own overlap.
- **RHF**: Löwdin orthogonalisation, LAPACKE `dsyev`, BLAS `dgemm`, Pulay DIIS.
- **Active space**: MO window `[first, first+norb)`; doubly occupied orbitals below the
  window fold into E_core and an effective one-body operator. Only closed-shell active
  spaces are supported.
- **Hamiltonian**: second-quantised H in spin orbitals (qubit `2p` = α, `2p+1` = β),
  mapped with Jordan–Wigner. Pauli strings are stored as `(x,z)` bitmasks.
- **Validation** (`fermion_validate`): integral symmetry, Hermitian Pauli coefficients,
  `‖[H,N]ψ‖` on random states, Hermiticity via random vectors, HF expectation equal to
  RHF.
- **FCI**: Fock-space determinant Hamiltonian with explicit sign rules (independent of
  the Jordan–Wigner mapping), diagonalised with `dsyev`. The JW Pauli Hamiltonian
  restricted to the N = 2 sector must reproduce it.
- **UCCSD**: spin-conserving singles and doubles (10 + 25 = 35 for CAS(2,6)). Each
  generator `T − T†` is written as Pauli strings and applied as first-order Trotter
  rotations `exp(iθ d P)`. This is exact at θ = 0 and exact for commuting strings; for
  this system it reaches FCI, but the Trotter error is not separately bounded.
- **Optimiser**: NLopt `NLOPT_LN_COBYLA` from θ = 0.
- **Circuits**: logical IR → OpenQASM 3 (`h, sdg, s, x, rz, cx`). Native compilation to
  `rz, sx, x, cz` (CX = H·CZ·H, H = RZ(π/2)·SX·RZ(π/2)). Routing inserts SWAPs
  (3 CX) along shortest paths inside the chosen layout.
- **Backend selection**: online QPUs with ≥ 12 qubits and `cz` in their basis; the
  candidate with the lowest layout error score (readout errors + CZ errors on the
  selected 12-qubit connected subgraph) is used unless `--backend` is given.
- **Hardware job**: Estimator V2, `program_id = "estimator"`, `version = 2`,
  `resilience_level = 1`, `default_precision` (default 0.02), observables as Pauli
  strings over the physical register (Qiskit label order, qubit 0 rightmost). The
  energy is Σ cᵢ⟨Pᵢ⟩ with σ = √(Σ cᵢ² σᵢ²), which ignores covariances between terms.

## Verification

`make test` runs:

- `test_chem`: basis normalisation and symmetry, ERI permutational symmetry, RHF and
  CAS(2,6) FCI against PySCF reference values, Hamiltonian validation, JW vs
  determinant FCI.
- `test_circuit`: Pauli rotations vs gate IR vs native decomposition; UCCSD counts,
  θ = 0 gives HF, 20 random states ≥ FCI, gate circuit == direct preparation,
  native circuit == IR.
- `test_ibm`: calibration parsing, layout selection avoids a noisy qubit, every native
  CZ on a coupler, and the routed 12-qubit state equals the logical state up to the
  final permutation. These tests use **synthetic in-memory devices**, not IBM data.

The PySCF numbers come from a scratch environment outside this repository. They are
hard-coded as reference values only.

## Status and limitations

- **No hardware results.** No `IBM_QUANTUM_API_KEY` or CRN was available, so no QPU job
  was submitted, and no hardware energy, uncertainty, or job ID exists. The REST calls
  follow the 2026-04-15 OpenAPI schema and the documented Estimator example, but they
  have not been exercised against the live service.
- **Routing is poor for this ansatz.** Jordan–Wigner strings of the UCCSD generators
  span many qubits. On a 12-qubit ring the routed circuit needs 4156 SWAPs; on a 4×4
  grid it needs 6752 SWAPs and about 22,600 native CZ gates. On current devices this
  depth will very likely destroy the signal. `hardware preflight` reports this before
  any submission.
- **Only CZ-basis backends** are supported. ECR-only devices are rejected.
- **Statevector limit**: 12 qubits, so only CAS(2,k) with k ≤ 6 can be validated
  locally.
- **Two-electron active spaces only for the UCCSD path.** FCI and the Hamiltonian are
  general in the electron count, but this build only exercises CAS(2,k).
- **Scope**: H₂ at a single geometry. No dissociation curve, no noise model, no
  error mitigation beyond the Estimator's `resilience_level = 1`.

## Layout

```
include/qchem.h          shared interfaces
src/common/              utilities, JSON results and metadata
src/basis/               molecule, NWChem basis parser, libcint environment
src/integrals/           S, T, V, ERI via libcint
src/scf/rhf.c            RHF (LAPACK + BLAS + DIIS)
src/active/active.c      AO→MO transform, active space, E_core
src/fermion/             Pauli algebra and Jordan–Wigner Hamiltonian
src/sim/statevec.c       12-qubit statevector simulator
src/qasm/circuit.c       circuit IR, Pauli compilation, native gates, OpenQASM 3
src/vqe/                 UCCSD ansatz and COBYLA driver
src/fci/fci.c            determinant FCI reference
src/hw/                  IBM REST client, calibration, layout/routing, hardware CLI
src/main.c               CLI
tests/                   scientific validation programs
scripts/fetch_libcint.sh libcint from conda-forge
data/basis/cc-pVQZ.nw    basis set (Basis Set Exchange)
```
