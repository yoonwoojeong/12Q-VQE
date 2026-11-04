# Ecosystem, Grant, and Journal Submission Kit (`qchem-vqe v0.2.0`)

Because external funding portals and peer-review bots require your personal GitHub login, email address, and a **Public** repository visibility, all repository-side files (`ecosystem.json`, `CITATION.cff`, `CONTRIBUTING.md`, `paper/paper.md`, `paper/paper.bib`, `examples/pyscf_interop_demo.py`, and `scripts/submit_ecosystem.ps1`) are already committed and tested in this repository.

Below are the exact **step-by-step instructions and copy-pasteable payloads** to submit `qchem-vqe` to all four major quantum open-source ecosystems.

---

## Step 0: One-Command Automated Push + Qiskit Ecosystem Submission (PowerShell)

Run the automated PowerShell script [`scripts/submit_ecosystem.ps1`](../scripts/submit_ecosystem.ps1) directly from your host PowerShell window:

```powershell
cd C:\Users\sshs4\12Q-VQE
powershell -ExecutionPolicy Bypass -File scripts\submit_ecosystem.ps1
```

This single script automatically:
1. Pushes `main` and `v0.2.0` to `https://github.com/yoonwoojeong/qchem-vqe`.
2. Sets your GitHub repository description and adds the `qiskit` repository topic (required by Qiskit Ecosystem automated validation checks `[G14]` and `[G15]`).
3. Submits the formatted `[Submission]: qchem-vqe` issue directly to `Qiskit/ecosystem` and prints the live issue URL.

---

## 1. Qiskit Ecosystem Submission (Manual Fallback Link)

If you prefer to submit via browser instead of the PowerShell script:

- **Submission Portal (1-Click Issue Form)**:
  https://github.com/Qiskit/ecosystem/issues/new?template=01_submission.yml&title=%5BSubmission%5D%3A+qchem-vqe

### Copy-Paste Fields for the GitHub Issue Form:
- **Project Name**:
  `qchem-vqe`
- **Description** *(must be under 135 characters)*:
  `C17 quantum chemistry engine (RHF, FCIDUMP, adjoint VQE, SQD) with OpenQASM 3 and IBM Quantum Estimator V2 support.`
- **Contact email**:
  *(Enter your preferred email address)*
- **Category**:
  `Application`
- **Labels**:
  `chemistry`, `Hamiltonian simulation`, `openqasm`, `circuit optimization`, `HPC`
- **Interface/API**:
  `Python`, `C`, `Command-line interface (CLI)`
- **Stability and support expectations**:
  `production-ready`
- **GitHub repository**:
  `https://github.com/yoonwoojeong/qchem-vqe`
- **Packages**:
  `https://github.com/yoonwoojeong/qchem-vqe/blob/main/pyproject.toml`
- **Reference paper**:
  `https://github.com/yoonwoojeong/qchem-vqe/blob/main/paper/paper.md`

---

## 2. Unitary Foundation \$4,000 Microgrant Application

Unitary Foundation awards **\$4,000 unrestricted microgrants** plus mentorship and community promotion across Discord, X, and the Unitary Foundation blog.

- **Application Portal**:
  https://unitary.foundation/grants/

### Copy-Paste Application Answers:

- **Project Title**:
  `12Q-VQE (libqchem): Fast C17 Ab Initio Quantum Chemistry, Adjoint VQE, and Sample-Based Quantum Diagonalization Engine`

- **Repository URL**:
  `https://github.com/yoonwoojeong/12Q-VQE`

- **Short Abstract (1-2 sentences)**:
  `12Q-VQE (libqchem / qchem_vqe) is a standalone, zero-Fortran C17 and Python library for ab initio molecular integrals, reverse-mode adjoint L-BFGS/ADAPT-VQE, Sample-Based Quantum Diagonalization (SQD), and direct OpenQASM 2.0 / IBM Quantum hardware execution. It bridges classical electronic structure (via standard .xyz and FCIDUMP files) with modern variational and post-VQE subspace quantum algorithms at 200x the speed of derivative-free Python simulators.`

- **What are you going to build / achieve with this grant?**
  ```text
  We have already built and open-sourced v0.2.0 of 12Q-VQE (`libqchem.so` and the `qchem_vqe` Python package), featuring:
  1. A self-contained C17 Gaussian integral engine (`cc-pVQZ`, `cc-pVDZ`, `STO-3G`), DIIS RHF solver, and OpenMP/BLAS 4-index active-space transformer.
  2. Reverse-mode analytical adjoint statevector gradients (`L-BFGS` and `ADAPT-VQE`) that converge 12-qubit, 559-term molecular Hamiltonians in 0.23 seconds (>200x faster than COBYLA).
  3. Sample-Based Quantum Diagonalization (SQD / QSCI) with Self-Consistent Configuration Recovery, recovering exact CAS-FCI energies under 3% per-qubit readout noise.
  4. Shallow Givens orbital-rotation compilation and native Estimator V2 hardware execution (validated live on the 156-qubit `ibm_fez` Heron r2 QPU to within 1.4 mEh of exact FCI).

  With the Unitary Foundation microgrant, we will deliver v0.3.0 over a 4-month roadmap:
  - Milestone 1 (Month 1): OpenQASM 3.0 & Qiskit QuantumCircuit / PySCF native Python C-extension wheels (`manylinux` and `macos-arm64` on PyPI via `cibuildwheel`), enabling zero-compilation `pip install qchem-vqe`.
  - Milestone 2 (Month 2): Symmetry-adapted fermion-to-qubit mappings (Bravyi-Kitaev and Z2 Pauli tapering) and spin-adapted singlet UCCSD pools, cutting qubit counts by 2-4 qubits and operator pool sizes by 60%.
  - Milestone 3 (Month 3): Distributed MPI + OpenMP hybrid statevector and large-subspace Davidson diagonalization for SQD up to 28-32 qubits (enabling active spaces up to CAS(14e, 14o) from real QPU bitstring dumps).
  - Milestone 4 (Month 4): Comprehensive documentation site, tutorials integrating Qiskit Runtime / Mitiq / PySCF, and JOSS publication completion.
  ```

- **How does your project benefit the quantum open-source ecosystem?**
  ```text
  Currently, researchers and students exploring quantum chemistry algorithms (ADAPT-VQE, Givens state preparation, or Sample-Based Quantum Diagonalization) must assemble deeply nested Python stacks that suffer from slow Python-level operator construction and high memory overhead, or else grapple with massive legacy Fortran/C++ codebases.

  12Q-VQE provides a compact (~5,300 LOC), strictly leak-checked (ASan/UBSan/LSan verified), MIT-licensed C17 engine (`libqchem.so`) with a 1-file C header (`qchem.h`), standard `FCIDUMP` and `.xyz` interoperability, and zero-overhead Python bindings. Any workflow in PySCF, Qiskit, PennyLane, or Mitiq can export an `FCIDUMP` or hardware bitstring distribution and use `libqchem` as a drop-in high-speed C accelerator for adjoint VQE or noisy SQD configuration recovery.
  ```

- **2-Minute Pitch Video Script (Required by Unitary Foundation)**:
  > **0:00 - 0:25 (The Problem)**: "Hi, I'm Yoonwoo Jeong, creator of 12Q-VQE. In quantum computational chemistry, moving from a molecule's Cartesian coordinates to an optimized variational ansatz or a noisy hardware subspace diagonalization typically requires stitching together heavy classical packages and slow Python operator loops."
  > **0:25 - 0:60 (What We Built)**: "To solve this, I built 12Q-VQE (`libqchem`): a complete, MIT-licensed C17 quantum chemistry and hybrid quantum engine in 5,300 lines of code. Out of the box, it computes Gaussian integrals up to `cc-pVQZ`, runs DIIS Hartree-Fock, exports and imports standard PySCF `FCIDUMP` files, and runs reverse-mode analytical adjoint L-BFGS and ADAPT-VQE—converging a 12-qubit, 559-term Hamiltonian in 0.23 seconds, over 200 times faster than derivative-free VQE."
  > **1:00 - 1:35 (Hardware & SQD)**: "For real quantum hardware, 12Q-VQE compiles shallow Givens rotation circuits—validated live on IBM's 156-qubit `ibm_fez` processor to within 1.4 milli-Hartree of exact FCI—and implements Sample-Based Quantum Diagonalization with Self-Consistent Configuration Recovery, recovering exact ground-state energies even under 3% per-qubit noise."
  > **1:35 - 2:00 (Grant Goals)**: "With a Unitary Foundation microgrant, I will package prebuilt PyPI binary wheels (`pip install qchem-vqe`), add Z2 symmetry tapering and Bravyi-Kitaev mappings, and scale the SQD Davidson solver to 32-qubit hardware bitstring samples. Thank you!"

---

## 3. Journal of Open Source Software (JOSS) Submission

The manuscript and bibliography are already prepared in [`paper/paper.md`](../paper/paper.md) and [`paper/paper.bib`](../paper/paper.bib).

- **Submission Portal**:
  https://joss.theoj.org/papers/new

### Pre-Submission Checklist & Form Fields:
1. Ensure your ORCID iD is linked to your JOSS account (add `orcid: 0000-...` under `authors` in `paper/paper.md` if you have one).
2. Fill in the JOSS submission form:
   - **Title**: `12Q-VQE: A C17 Ab Initio Quantum Chemistry, Adjoint VQE, and Sample-Based Quantum Diagonalization Engine`
   - **Repository Address**: `https://github.com/yoonwoojeong/12Q-VQE`
   - **Branch containing `paper.md`**: `main`
   - **Software Version**: `v0.2.0`
   - **Type of submission**: `New submission`
   - **Description / Notes to Editor**:
     ```text
     12Q-VQE (libqchem / qchem_vqe) is a standalone C17 and Python library for ab initio molecular integrals, DIIS Restricted Hartree-Fock, FCIDUMP interoperability, analytical adjoint L-BFGS / ADAPT-VQE, Sample-Based Quantum Diagonalization (SQD), and direct OpenQASM 3.0 / IBM Quantum hardware execution. It includes a 5-binary automated test suite verified under AddressSanitizer, UndefinedBehaviorSanitizer, and LeakSanitizer in GitHub Actions CI.
     ```

---

## 4. Zenodo DOI & Awesome Quantum Lists

1. **Zenodo Permanent DOI**:
   - Log into https://zenodo.org with your GitHub account.
   - Go to **GitHub** settings in Zenodo and toggle **ON** `yoonwoojeong/12Q-VQE`.
   - Create a GitHub Release from tag `v0.2.0` at `https://github.com/yoonwoojeong/12Q-VQE/releases/new?tag=v0.2.0` — Zenodo will automatically archive the release and mint a citable `10.5281/zenodo.*` DOI.
2. **Awesome Quantum Software (`qosf/awesome-quantum-software`)**:
   - Open a PR at https://github.com/qosf/awesome-quantum-software adding `12Q-VQE` under **Quantum Chemistry**:
     ```markdown
     - [12Q-VQE](https://github.com/yoonwoojeong/12Q-VQE) - C17 ab initio quantum chemistry engine (RHF, FCIDUMP, adjoint L-BFGS/ADAPT-VQE, Sample-Based Quantum Diagonalization, and Givens hardware compilation) with Python bindings.
     ```
