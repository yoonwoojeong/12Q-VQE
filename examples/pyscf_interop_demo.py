#!/usr/bin/env python3
"""
examples/pyscf_interop_demo.py
==============================
End-to-end demonstration of `qchem_vqe` (`12Q-VQE`) as both:
1. A zero-dependency C17 quantum chemistry + VQE + SQD engine in Python
2. A drop-in `FCIDUMP` / `.xyz` accelerator for PySCF and Qiskit workflows

Run directly from the repository root:
    PYTHONPATH=python python3 examples/pyscf_interop_demo.py
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

# Ensure local python package is importable when running from repo root
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT / "python") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "python"))

import qchem_vqe as qc


def demo_xyz_and_fcidump_workflow() -> None:
    print("=" * 72)
    print("1. Arbitrary .xyz -> RHF -> CAS Active Space -> FCIDUMP Export/Import")
    print("=" * 72)

    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_path = os.path.join(tmpdir, "beh2_stretched.xyz")
        fcidump_path = os.path.join(tmpdir, "beh2_cas44.fcidump")

        # Write a standard Cartesian .xyz file (Angstroms)
        with open(xyz_path, "w", encoding="utf-8") as f:
            f.write("3\n")
            f.write("Linear BeH2 at R(Be-H) = 1.326 Angstrom\n")
            f.write("Be  0.000000  0.000000  0.000000\n")
            f.write("H   0.000000  0.000000  1.326000\n")
            f.write("H   0.000000  0.000000 -1.326000\n")

        # Run RHF + CAS(4e, 4o) in sto-3g and export standard Molpro/PySCF FCIDUMP
        cas = qc.export_fcidump(
            fcidump_path,
            xyz=xyz_path,
            basis="sto-3g",
            electrons=4,
            orbitals=4,
        )
        print(f"  Loaded XYZ     : {xyz_path}")
        print(f"  Exported       : {fcidump_path} ({os.path.getsize(fcidump_path)} bytes)")
        print(f"  Qubits (2*M)   : {cas['qubits']} ({cas['pauli_terms']} Pauli terms)")
        print(f"  RHF Energy     : {cas['rhf_full_basis_Eh']:.10f} Eh")
        print(f"  Exact CAS-FCI  : {cas['fci_cas_Eh']:.10f} Eh")

        # Re-import the FCIDUMP and solve with exact FCI, L-BFGS VQE, and SQD
        print("\n" + "=" * 72)
        print("2. Solving Imported FCIDUMP with Exact FCI, Adjoint L-BFGS VQE, and SQD")
        print("=" * 72)

        vqe_res = qc.run_vqe(
            fcidump=fcidump_path,
            ansatz="uccsd",
            optimizer="lbfgs",
        )
        sqd_res = qc.run_sqd(
            fcidump=fcidump_path,
            shots=4096,
            noise=0.02,
            sc_iters=4,
        )

        print(f"  Exact CAS-FCI  : {vqe_res['fci_cas_Eh']:.10f} Eh")
        print(
            f"  L-BFGS VQE     : {vqe_res['vqe_local_Eh']:.10f} Eh "
            f"(|dE| = {abs(vqe_res['vqe_minus_fci_Eh']) * 1e3:.4f} mEh, "
            f"{vqe_res['seconds'] * 1e3:.1f} ms)"
        )
        print(
            f"  SQD (2% noise) : {sqd_res['sqd_energy_Eh']:.10f} Eh "
            f"(|dE| = {abs(sqd_res['sqd_minus_fci_Eh']) * 1e3:.4f} mEh, "
            f"subspace = {sqd_res['subspace_dim']}/{sqd_res['fci_dim']} dets)"
        )


def demo_optional_pyscf_bridge() -> None:
    print("\n" + "=" * 72)
    print("3. Optional PySCF Interoperability Check")
    print("=" * 72)
    try:
        from pyscf import gto, scf  # type: ignore
        from pyscf.tools import fcidump  # type: ignore
    except ImportError:
        print(
            "  [Info] PySCF is not installed in this Python environment.\n"
            "  When PySCF is installed, any `pyscf.tools.fcidump.from_scf(mf, 'mol.fcidump')`\n"
            "  file can be passed directly to `qchem_vqe.run_vqe(fcidump='mol.fcidump')`\n"
            "  or `qchem_vqe.run_sqd(fcidump='mol.fcidump')`."
        )
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        fcidump_path = os.path.join(tmpdir, "pyscf_h2.fcidump")
        mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", verbose=0)
        mf = scf.RHF(mol).run()
        fcidump.from_scf(mf, fcidump_path)
        res = qc.run_vqe(fcidump=fcidump_path, ansatz="adapt", optimizer="lbfgs")
        print(f"  PySCF RHF (H2/sto-3g)  : {mf.e_tot:.10f} Eh")
        print(f"  12Q-VQE ADAPT-VQE      : {res['vqe_local_Eh']:.10f} Eh (FCI: {res['fci_cas_Eh']:.10f} Eh)")


if __name__ == "__main__":
    demo_xyz_and_fcidump_workflow()
    demo_optional_pyscf_bridge()
