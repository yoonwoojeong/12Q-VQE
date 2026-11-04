# Contributing to 12Q-VQE

Thank you for your interest in contributing to **12Q-VQE** (`libqchem` / `qchem_vqe`)! We welcome bug reports, new basis sets, quantum algorithm modules, hardware backend adapters, and documentation improvements.

---

## Development Setup

### Prerequisites
- A C17-compliant compiler (`gcc >= 10` or `clang >= 12`)
- OpenMP (`libomp-dev` / `libgomp`)
- BLAS (`libopenblas-dev` or Accelerate on macOS)
- `libcurl` (`libcurl4-openssl-dev`) for IBM Quantum Runtime REST commands
- Python 3.8+ (optional, for `qchem_vqe` bindings and interoperability tests)

### Building & Running the Test Suite

```bash
# 1. Build optimized CLI and shared library (build/qchem_vqe and build/libqchem.so)
make all

# 2. Run the 5-binary unit and integration test suite
make test

# 3. Run AddressSanitizer + UndefinedBehaviorSanitizer + LeakSanitizer verification
make sanitize

# 4. Run the Python ctypes interoperability suite
PYTHONPATH=python python3 examples/pyscf_interop_demo.py
```

---

## Code Standards & Architecture

1. **Zero Hidden Dependencies in Core C**:
   - The core library (`src/basis/`, `src/hf/`, `src/active/`, `src/mapping/`, `src/vqe/`, `src/fci/`, `src/common/`) must compile cleanly under `-std=c17 -Wall -Wextra -Werror`.
   - Use standard C17, OpenMP, and CBLAS primitives. Do not introduce heavy C++ or Fortran runtime dependencies into `libqchem`.
2. **Memory Safety**:
   - Every pull request must pass `make sanitize` (`-fsanitize=address,undefined` with LeakSanitizer enabled) with **zero leaks** and **zero undefined behavior**.
3. **Physical & Numerical Verification**:
   - New integrals, fermion-to-qubit mappings, or eigensolvers must include a standalone test in `tests/` validating results against exact analytical, FCI, or PySCF reference values (typically to `< 1e-8 Eh` for integrals/RHF and `< 1.6 mEh` chemical accuracy for VQE/SQD).
4. **C ABI & Python Bindings**:
   - Public C functions should be declared in `include/qchem.h` and automatically exported into `build/libqchem.so`. If a feature is user-facing, expose a clean wrapper in `python/qchem_vqe/__init__.py`.

---

## Pull Request Checklist

- [ ] `make all` compiles with zero warnings.
- [ ] `make test` passes all unit tests.
- [ ] `make sanitize` passes with zero ASan/UBSan/LSan errors.
- [ ] Added or updated unit tests in `tests/` covering the new functionality.
- [ ] Updated `README.md` or `include/qchem.h` docstrings if public APIs changed.
