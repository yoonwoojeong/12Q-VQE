# qchem-vqe (`libqchem`)

**Native C17 Quantum Chemistry VQE, ADAPT-VQE, and Sample-Based Quantum Diagonalization (SQD) Engine with Direct IBM Quantum QPU Integration.**

`qchem-vqe` is a zero-Python, high-performance electronic structure and hybrid quantum-classical engine written in portable C17. It computes ab initio Gaussian atomic orbital integrals via [`libcint`](https://github.com/sunqm/libcint), solves full-basis Restricted Hartree-Fock (RHF) with Pulay DIIS, transforms active spaces (`CAS(N_e, N_o)`) with BLAS `dgemm`, maps second-quantized Hamiltonians via Jordan-Wigner, and solves active-space quantum chemistry via:
- **Exact Determinant CAS-FCI** (OpenMP-parallelized Slater-Condon subspace solver)
- **UCCSD & ADAPT-VQE** accelerated by **exact reverse-mode adjoint gradients** + **L-BFGS** (`>200x` faster than derivative-free COBYLA)
- **Sample-Based Quantum Diagonalization (SQD / QSCI)** with **Self-Consistent Configuration Recovery (S-CORE)** for noise-resilient chemical accuracy on noisy quantum hardware
- **Shallow Natural-Orbital Givens Compilation & 1D Path Routing** (`0 SWAPs`, native `rz, sx, x, cz` gates in OpenQASM 3.0)
- **Direct IBM Quantum Compute Service REST Client** (`IBM-API-Version: 2026-04-15`, Estimator V2, multi-PUB controlled experiments, live calibration snapshots)
- **Standard Ecosystem Interoperability** via **`FCIDUMP`** (PySCF, ORCA, Molpro, Psi4) and **`.xyz`** molecular geometry files, plus a zero-dependency Python wrapper (`qchem_vqe`).

---

## Key Capabilities & Benchmarks

| System | Basis / Active Space | Qubits | Method | Circuit 2Q Gates | Wall Time | Error vs CAS-FCI ($\text{E}_\text{h}$) |
| :--- | :--- | :---: | :--- | :---: | :---: | :---: |
| $\text{H}_2$ ($0.735\text{ \AA}$) | `cc-pVQZ` `CAS(2,6)` | 12 | UCCSD + COBYLA | 2,400 `CX` | `48.80 s` | $6.98 \times 10^{-11}$ |
| $\text{H}_2$ ($0.735\text{ \AA}$) | `cc-pVQZ` `CAS(2,6)` | 12 | **UCCSD + L-BFGS (Adjoint Grad)** | 2,400 `CX` | **`0.23 s`** | **$1.33 \times 10^{-15}$** |
| $\text{H}_2$ ($0.735\text{ \AA}$) | `cc-pVQZ` `CAS(2,6)` | 12 | **ADAPT-VQE (9 ops) + L-BFGS** | **592 `CX`** | **`0.66 s`** | **$2.44 \times 10^{-15}$** |
| $\text{H}_2$ ($0.735\text{ \AA}$) | `cc-pVQZ` `CAS(2,6)` | 12 | **Givens Natural-Orbital (`r=5`)** | **368 native `CZ` (`0 SWAP`)** | **`0.01 s`** | **$3.96 \times 10^{-6}$** |
| $\text{H}_2$ ($0.735\text{ \AA}$) | `cc-pVQZ` `CAS(2,8)` | 16 | **16-Qubit CAS-FCI (1,241 Paulis)** | — | **`0.08 s`** | Exact ($-1.14944547\text{ E}_\text{h}$) |
| $\text{LiH}$ ($1.595\text{ \AA}$) | `sto-3g` `CAS(2,5)` (frozen $1s$) | 10 | **ADAPT-VQE (10 ops) + L-BFGS** | **528 `CX`** | **`0.89 s`** | **$2.66 \times 10^{-14}$** |
| $\text{H}_4$ ($1.000\text{ \AA}$) | `sto-3g` `CAS(4,4)` ($4e^-$) | 8 | **ADAPT-VQE (8 ops) + L-BFGS** | **400 `CX`** | **`0.40 s`** | **$3.10 \times 10^{-10}$** |
| $\text{BeH}_2$ ($1.326\text{ \AA}$) | `sto-3g` `CAS(4,6)` ($4e^-$) | 12 | **ADAPT-VQE (18 / 92 ops)** | **992 `CX`** | **`2.92 s`** | **$2.68 \times 10^{-6}$** |

### Noise Resilience: SQD / QSCI vs Raw Estimator
Under `3.0%` per-qubit depolarizing/readout bit-flip noise (where `>25%` of all shots leak out of the particle-number sector and raw Pauli `Estimator` expectations lose `>0.16 Eh`):
- **$\text{H}_2$ `CAS(2,6)` (12 qubits, 4,096 shots)**: Raw noisy Estimator error is **$+0.165198\text{ E}_\text{h}$**; **SQD with S-CORE** recovers **$-1.148902977328\text{ E}_\text{h}$** (error **$2.22 \times 10^{-16}\text{ E}_\text{h}$**, using 22 / 66 determinants in `94 ms`).
- **$\text{LiH}$ `CAS(2,5)` (`FCIDUMP`, 10 qubits, 4,096 shots)**: Raw noisy Estimator error is **$+0.085043\text{ E}_\text{h}$**; **SQD with S-CORE** recovers **$-7.554736489187\text{ E}_\text{h}$** (error **$8.88 \times 10^{-16}\text{ E}_\text{h}$**, using 16 / 45 determinants in `63 ms`).

---

## Installation & Deployment

### Option 1: Build from Source (Linux / WSL Ubuntu)
Requirements: GCC or Clang (C17 + OpenMP), `liblapacke-dev`, `libopenblas-dev`, `libnlopt-dev`, `libcjson-dev`, `libcurl4-openssl-dev`, `zstd`, `unzip`, `curl`.

```bash
make fetch-libcint          # downloads libcint 6.1.3 from conda-forge into third_party/libcint
make -j4                    # builds build/qchem-vqe, build/libqchem.a, and build/libqchem.so
make test                   # runs all 5 scientific & hardware validation suites
make sanitize               # runs full suite under AddressSanitizer + UndefinedBehaviorSanitizer + LeakSanitizer
sudo make install           # installs to /usr/local (bin, lib, include, share/qchem-vqe, pkg-config qchem.pc)
```

### Option 2: Portable Binary Release Tarball
Download `qchem-vqe-linux-x86_64.tar.gz` from [GitHub Releases](https://github.com/yoonwoojeong/12Q-VQE/releases) (or build it locally with `make dist`):
```bash
tar -xzf qchem-vqe-linux-x86_64.tar.gz
./qchem-vqe/bin/qchem-vqe sqd --mol H2 --basis cc-pVQZ --shots 4096
```

### Option 3: Docker / Apptainer Container
```bash
docker build -t qchem-vqe .
docker run --rm qchem-vqe sqd --mol LiH --basis sto-3g --electrons 2 --orbitals 5
```

### Option 4: Python Bridge (`pip install`)
```bash
pip install -e .
python3 -c "import qchem_vqe; print(qchem_vqe.run_sqd(mol='LiH', basis='sto-3g', electrons=2, orbitals=5))"
```

### Option 5: Windows (PowerShell via WSL)
```powershell
powershell -ExecutionPolicy Bypass -File scripts\quickstart-wsl.ps1
```

---

## Command-Line Interface

```bash
# 1. Integral preparation & full-basis RHF (built-in molecules: H2, HeH+, LiH, BeH2, H4, or custom --xyz FILE)
./build/qchem-vqe prepare --mol LiH --basis sto-3g --bond 1.595
./build/qchem-vqe rhf --mol H2 --basis cc-pVQZ --bond 0.735
./build/qchem-vqe rhf --xyz my_molecule.xyz --charge 0 --mult 1 --basis cc-pVDZ

# 2. Active-Space Hamiltonian, Validation, FCI, and FCIDUMP Export/Import
./build/qchem-vqe active-space --mol H2 --basis cc-pVQZ --electrons 2 --orbitals 6
./build/qchem-vqe active-space --mol LiH --basis sto-3g --electrons 2 --orbitals 5 --fcidump-out lih.fcidump
./build/qchem-vqe active-space --fcidump lih.fcidump

# 3. Local VQE: L-BFGS with Reverse-Mode Adjoint Gradients (default) or ADAPT-VQE
./build/qchem-vqe vqe --mode local --mol H2 --basis cc-pVQZ --electrons 2 --orbitals 6
./build/qchem-vqe vqe --mode local --mol BeH2 --basis sto-3g --electrons 4 --orbitals 6 --ansatz adapt

# 4. Sample-Based Quantum Diagonalization (SQD / QSCI) with Configuration Recovery
./build/qchem-vqe sqd --mol H2 --basis cc-pVQZ --electrons 2 --orbitals 6 --shots 4096 --noise 0.03
./build/qchem-vqe sqd --fcidump lih.fcidump --shots 4096 --noise 0.03 --sc-iters 4

# 5. Potential Energy Surface (PES) Dissociation Curve Scans
./build/qchem-vqe pes --mol H2 --basis cc-pVDZ --electrons 2 --orbitals 4 --rmin 0.5 --rmax 2.5 --npts 9
./build/qchem-vqe pes --mol LiH --basis sto-3g --electrons 2 --orbitals 5 --rmin 1.0 --rmax 2.8 --npts 7

# 6. IBM Quantum Hardware Inspection, Preflight, Submission & Controlled Experiments
./build/qchem-vqe auth login                                  # hidden prompt for IBM Cloud API key & CRN
./build/qchem-vqe hardware list
./build/qchem-vqe hardware preflight --backend ibm_fez
./build/qchem-vqe vqe --mode hardware --backend ibm_fez --submit
./build/qchem-vqe job status JOB_ID
./build/qchem-vqe job results JOB_ID
./build/qchem-vqe experiment plan experiments/h2_cas26.ini
./build/qchem-vqe experiment run experiments/h2_cas26.ini
./build/qchem-vqe experiment results JOB_ID
```

---

## Live IBM Quantum Hardware Record (`ibm_fez`)

A live 12-qubit Estimator V2 job (`dblcn1pkbkrc73cnfv4g`) was executed on the 156-qubit **`ibm_fez`** Heron processor using the truncated natural-orbital Givens circuit (`r = 5` pairs, 327 Pauli observables):
- **Physical 1D path**: `[1, 2, 3, 16, 23, 24, 25, 37, 45, 44, 43, 56]`
- **Compiled circuit**: `1,250` native gates, **`368` native `CZ` gates**, **`0` routing SWAPs**, critical-path duration **`30.0 us`** (vs worst-case $T_2 = 61\text{ }\mu\text{s}$ on the layout)
- **Artifacts**: [`results/hardware/dblcn1pkbkrc73cnfv4g.json`](results/hardware/dblcn1pkbkrc73cnfv4g.json), [`results/hardware/dblcn1pkbkrc73cnfv4g.qasm`](results/hardware/dblcn1pkbkrc73cnfv4g.qasm), and [`results/hardware/dblcn1pkbkrc73cnfv4g_results.json`](results/hardware/dblcn1pkbkrc73cnfv4g_results.json).
- Because unmitigated 368-`CZ` expectation values suffer exponential Pauli damping across 327 terms, this hardware benchmark directly motivates the **controlled angle-zero noise subtraction (`experiment`)** and **Sample-Based Quantum Diagonalization (`sqd`)** workflows built into `qchem-vqe`.

---

## Architecture & Repository Layout

```text
include/qchem.h            Public C17 library API (libqchem)
src/common/                Memory-safe helpers, dynamic path resolution, JSON metadata
src/basis/                 Multi-element molecules, .xyz loader, NWChem S/P/D/F/G + SP parser
src/integrals/             Spherical 1e/2e AO integrals via libcint
src/scf/rhf.c              Full-basis RHF solver (Löwdin orthogonalization + Pulay DIIS)
src/active/active.c        BLAS dgemm AO->MO integral transform, frozen core, FCIDUMP I/O
src/fermion/               Bitmask Pauli algebra and Jordan-Wigner Hamiltonian + validation
src/sim/statevec.c         OpenMP 20-qubit statevector simulator and bit-twiddling Pauli kernels
src/qasm/circuit.c         Logical IR, chain CNOT ladders, native rz/sx/x/cz compiler, OpenQASM 3.0
src/vqe/uccsd.c            Spin-conserving UCCSD excitation pool
src/vqe/vqe.c              Reverse-mode adjoint gradients, L-BFGS / COBYLA, and ADAPT-VQE
src/vqe/givens.c           Two-electron SVD natural-orbital Givens circuit compiler
src/fci/fci.c              Subspace Slater-Condon FCI solver and SQD / S-CORE recovery engine
src/hw/                    IBM Quantum REST client, 1D path router, hardware CLI, controlled experiments
python/qchem_vqe/          Zero-dependency Python bridge (ctypes + CLI JSON wrapper)
tests/                     5 automated test suites (chem, circuit, creds, experiment, ibm)
data/basis/                NWChem basis sets (STO-3G, cc-pVDZ, cc-pVQZ)
examples/                  End-to-end Python & FCIDUMP/PySCF interoperability scripts
paper/                     Journal of Open Source Software (JOSS) manuscript (paper.md, paper.bib)
docs/                      Ecosystem, grant, and journal submission guides
```

---

## Citation & Community

- **Cite this software**: See [`CITATION.cff`](CITATION.cff) or the JOSS manuscript in [`paper/paper.md`](paper/paper.md).
- **Contributing**: See [`CONTRIBUTING.md`](CONTRIBUTING.md) for code standards and sanitizer requirements.
- **Ecosystem & Grant Submissions**: See [`docs/ECOSYSTEM_SUBMISSIONS.md`](docs/ECOSYSTEM_SUBMISSIONS.md) for Qiskit Ecosystem, Unitary Foundation, and JOSS submission details.
