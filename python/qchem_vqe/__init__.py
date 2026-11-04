"""
qchem_vqe: Zero-dependency Python bridge to the native C17 libqchem / qchem-vqe engine.

Supports:
  - Direct ctypes binding to libqchem.so (FCIDUMP / FCI / SQD / VQE)
  - High-level CLI/JSON wrappers for VQE, ADAPT-VQE, PES scans, and SQD (Sample-Based Quantum Diagonalization)
  - Seamless interoperability with PySCF and Qiskit via FCIDUMP and .xyz files.
"""

from __future__ import annotations
import ctypes
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

__version__ = "0.2.0"


def _find_repo_root() -> Optional[Path]:
    env_root = os.environ.get("QCHEM_ROOT_DIR")
    if env_root and Path(env_root).is_dir():
        return Path(env_root)
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "build" / "qchem-vqe").is_file() or (parent / "include" / "qchem.h").is_file():
            return parent
    return None


def find_binary() -> str:
    root = _find_repo_root()
    if root and (root / "build" / "qchem-vqe").is_file():
        return str(root / "build" / "qchem-vqe")
    which = shutil.which("qchem-vqe")
    if which:
        return which
    raise FileNotFoundError(
        "qchem-vqe binary not found. Run `make all` in the repository root or `sudo make install`."
    )


def find_library() -> ctypes.CDLL:
    root = _find_repo_root()
    candidates = []
    if root:
        candidates.append(root / "build" / "libqchem.so")
    candidates.extend([
        Path("/usr/local/lib/libqchem.so"),
        Path("/usr/lib/libqchem.so"),
    ])
    for cand in candidates:
        if cand.is_file():
            return ctypes.CDLL(str(cand))
    raise FileNotFoundError("libqchem.so not found. Run `make all` or `sudo make install`.")


def run_vqe(
    mol: str = "H2",
    basis: str = "cc-pVQZ",
    electrons: int = 2,
    orbitals: int = 6,
    bond: float = 0.735,
    ansatz: str = "adapt",
    optimizer: str = "lbfgs",
    fcidump: Optional[str] = None,
    xyz: Optional[str] = None,
) -> Dict[str, Any]:
    """Run local VQE (UCCSD or ADAPT-VQE with analytic adjoint gradients) and return parsed JSON results."""
    bin_path = find_binary()
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            bin_path, "vqe", "--mode", "local",
            "--mol", mol,
            "--basis", basis,
            "--electrons", str(electrons),
            "--orbitals", str(orbitals),
            "--bond", str(bond),
            "--ansatz", ansatz,
            "--opt", optimizer,
            "--results", tmpdir,
        ]
        if fcidump:
            cmd.extend(["--fcidump", str(fcidump)])
        if xyz:
            cmd.extend(["--xyz", str(xyz)])
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        json_files = list(Path(tmpdir).glob("vqe_local_CAS*.json"))
        if not json_files:
            raise RuntimeError("qchem-vqe did not produce a VQE JSON result file.")
        return json.loads(json_files[0].read_text())


def run_sqd(
    mol: str = "H2",
    basis: str = "cc-pVQZ",
    electrons: int = 2,
    orbitals: int = 6,
    bond: float = 0.735,
    shots: int = 4096,
    noise: float = 0.02,
    sc_iters: int = 4,
    fcidump: Optional[str] = None,
) -> Dict[str, Any]:
    """Run Sample-Based Quantum Diagonalization (SQD / QSCI) with self-consistent configuration recovery."""
    bin_path = find_binary()
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            bin_path, "sqd",
            "--mol", mol,
            "--basis", basis,
            "--electrons", str(electrons),
            "--orbitals", str(orbitals),
            "--bond", str(bond),
            "--shots", str(shots),
            "--noise", str(noise),
            "--sc-iters", str(sc_iters),
            "--results", tmpdir,
        ]
        if fcidump:
            cmd.extend(["--fcidump", str(fcidump)])
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        json_files = list(Path(tmpdir).glob("sqd_*.json"))
        if not json_files:
            raise RuntimeError("qchem-vqe did not produce an SQD JSON result file.")
        return json.loads(json_files[0].read_text())


def run_pes(
    mol: str = "H2",
    basis: str = "cc-pVDZ",
    electrons: int = 2,
    orbitals: int = 4,
    rmin: float = 0.5,
    rmax: float = 2.5,
    npts: int = 9,
) -> Dict[str, Any]:
    """Run a Potential Energy Surface (PES) dissociation scan and return all RHF, FCI, UCCSD, and ADAPT-VQE points."""
    bin_path = find_binary()
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            bin_path, "pes",
            "--mol", mol,
            "--basis", basis,
            "--electrons", str(electrons),
            "--orbitals", str(orbitals),
            "--rmin", str(rmin),
            "--rmax", str(rmax),
            "--npts", str(npts),
            "--results", tmpdir,
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        json_files = list(Path(tmpdir).glob("pes_*.json"))
        if not json_files:
            raise RuntimeError("qchem-vqe did not produce a PES JSON result file.")
        return json.loads(json_files[0].read_text())
