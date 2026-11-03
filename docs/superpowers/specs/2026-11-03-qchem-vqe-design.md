# qchem-vqe: H2 / CAS(2,6) VQE in C17 with IBM Quantum REST — design

Status: classified **architectural**. The request explicitly asks for autonomous
implementation with decisions made by the implementer, so the spec is written
and implemented in the same turn; decisions below are listed so they can be
redirected.

## Goal

A C17 command-line program that
1. computes cc-pVQZ AO integrals for H2 (0.735 Å) with libcint,
2. runs RHF in the full cc-pVQZ basis,
3. builds a configurable active space (default CAS(2,6): 2 electrons, 6 spatial
   orbitals = 12 spin orbitals = 12 qubits) via AO→MO transformation,
4. maps the active-space Hamiltonian to Paulis with Jordan–Wigner, validating
   Hermiticity and particle-number conservation,
5. runs spin-conserving UCCSD (Trotterised Pauli rotations) on a 12-qubit
   statevector simulator with COBYLA (NLopt),
6. computes a determinant-based CAS-FCI reference,
7. talks to the IBM Quantum Compute Service REST API (IAM, backends,
   calibration, Estimator V2 jobs) and runs the optimised circuit on a real QPU
   only after an explicit interactive confirmation.

## Verified external facts

- IBM REST API version used: `IBM-API-Version: 2026-04-15` (latest in the docs).
  Base URL `https://quantum.cloud.ibm.com/api/v1` (EU: `https://eu-de.quantum.cloud.ibm.com/api/v1`).
- Auth: `POST https://iam.cloud.ibm.com/identity/token`,
  `application/x-www-form-urlencoded`, `grant_type=urn:ibm:params:oauth:grant-type:apikey&apikey=…`.
  Headers on every call: `Authorization: Bearer <token>`, `Service-CRN: <instance CRN>`, `IBM-API-Version`.
- Endpoints used: `GET /backends`, `GET /backends/{id}/configuration|properties|status`,
  `POST /jobs`, `GET /jobs/{id}`, `GET /jobs/{id}/results`, `POST /jobs/{id}/cancel`.
- Estimator job body: `{"program_id":"estimator","backend":NAME,"params":{"pubs":[[qasm3, [pauli labels…]]],"version":2,"resilience_level":1,"options":{…}}}`.
  Observables are passed as Pauli label strings (rightmost character = qubit 0).
  Results: `data[i].results.{evs,stds}` as `CompressedTensorModel` (base64 little-endian f64, C order).
- libcint 6.1.3 (conda-forge build; distro packages do not provide libcint).
  Env layout: `PTR_EXPCUTOFF=0`, `PTR_COMMON_ORIG=1`, `PTR_ENV_START=20`,
  atm slots `CHARGE_OF=0, PTR_COORD=1, NUC_MOD_OF=2`,
  bas slots `ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF`.
- Basis: cc-pVQZ for H from Basis Set Exchange (NWChem format), `[4s3p2d1f]` = 30
  spherical functions per H, 60 AO for H2.

## Decisions

| Topic | Decision |
|---|---|
| Qubit ordering | qubit `2p` = α spin-orbital of active MO `p`, qubit `2p+1` = β. Bit k of a basis index = occupation of qubit k. HF state = index 3. |
| Jordan–Wigner | `a_j = (∏_{k<j} Z_k)(X_j + iY_j)/2`. Pauli strings stored as `(x,z)` bitmasks with canonical operator `i^{popcount(x&z)} X^x Z^z` (Hermitian). |
| Active space | Window of `--orbitals` MOs starting at `first` (default: number of frozen-by-window core orbitals so that the window covers the active electrons). Inactive doubly occupied MOs fold into `E_core` and an effective 1-body operator. Only 2 active electrons are supported by the UCCSD/FCI path. |
| ERI storage | Full AO ERI tensor (60⁴ doubles ≈ 100 MB). Simple and fast enough for H2. |
| SCF | Löwdin orthogonalisation, LAPACKE `dsyev`, Pulay DIIS on the Fock matrix. Convergence: |ΔE| < 1e-10, density RMS < 1e-8. |
| UCCSD | Spin-conserving singles and doubles from the occupied spin orbitals. Generator `T−T†` is decomposed to Pauli strings; each string is a rotation `exp(iθ d_s P_s)`, first-order Trotter across the list (exact when the strings of an excitation commute). 35 parameters for CAS(2,6). |
| Optimiser | NLopt `NLOPT_LN_COBYLA`, start at θ=0 (HF), initial step 0.1, tolerance 1e-8. |
| FCI reference | Fock-space determinant Hamiltonian built from spin-orbital integrals via explicit creation/annihilation sign rules (independent of JW). Dense `dsyev`. |
| Cross-check | Lowest eigenvalue of the JW Pauli Hamiltonian restricted to the N=2 subspace must equal the determinant FCI. |
| Circuits | Logical IR → OpenQASM 3. Pauli rotations compile to H/S†-basis changes, CNOT ladders and RZ. Native compile for CZ-based backends: `rz, sx, x, cz` (H = RZ(π/2)·SX·RZ(π/2), CX = H·CZ·H). ECR-only backends are rejected. |
| Routing | Greedy: for each non-adjacent CX, SWAP the control along the shortest physical path (SWAP = 3 CX). Initial layout: connected 12-qubit subgraph chosen by summed CZ error + readout error. |
| Backend selection | Default: online backends with ≥12 qubits and `cz` in basis gates, ranked by the layout score. `--backend NAME` overrides. |
| Hardware submission | `vqe --mode hardware --submit` prints the plan (backend, layout, depth, gate count, #observables, estimated shots) and requires typing `SUBMIT`, unless `--yes` is given. |
| Credentials | `IBM_QUANTUM_API_KEY`, `IBM_QUANTUM_CRN`, optional `IBM_QUANTUM_API_URL`. Values are never printed or written to results. |
| Results | JSON files under `results/` containing numbers, versions, git revision, timestamps, and job IDs. |

## Validation plan

- Basis/integrals: overlap diagonal = 1, symmetry of S/T/V/ERI, comparison with
  an external PySCF run done in a scratch environment (not part of the build).
- RHF: energy compared with PySCF RHF/cc-pVQZ at the same geometry.
- Hamiltonian: Hermiticity, real Pauli coefficients, `[H,N]=0` on random vectors,
  HF expectation equals RHF active-space energy.
- FCI: determinant FCI vs Pauli-restricted diagonalisation.
- Simulator: unitarity/norm, gate-level IR vs direct Pauli rotations, decomposition
  to native gates equal up to global phase.
- VQE: energy at θ=0 equals HF; VQE energy ≥ FCI − 1e-10 (variational).
- IBM layer: parsing and routing are tested against small synthetic JSON fixtures
  clearly labelled as fixtures. No hardware numbers are fabricated; hardware runs
  happen only if credentials and explicit confirmation are provided.
