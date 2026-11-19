# Relativistic Particle Accelerator Simulator

[![Test, Build WASM & Deploy GitHub Pages](https://github.com/Mr-Kumar-Abhishek/Relativistic-Particle-Accelerator-Simulator/actions/workflows/deploy.yml/badge.svg)](https://github.com/Mr-Kumar-Abhishek/Relativistic-Particle-Accelerator-Simulator/actions/workflows/deploy.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

A high-accuracy, double-precision (`f64`) relativistic charged-particle accelerator simulator written in Rust (`particle-accelarator` crate) with zero external dependencies, validated via a 5-Sprint Test-Driven Development (TDD) physics suite, and compiled to WebAssembly (`wasm32-unknown-unknown`) for interactive browser visualization on GitHub Pages.

- **Repository:** https://github.com/Mr-Kumar-Abhishek/Relativistic-Particle-Accelerator-Simulator
- **Live WebAssembly Simulator:** https://mr-kumar-abhishek.github.io/Relativistic-Particle-Accelerator-Simulator/
- **License:** [BSD 3-Clause License](LICENSE)

---

## Physics & Numerical Accuracy Highlights
- **CODATA 2022 Exact Physical Constants (`f64`):** Speed of light $c$, elementary charge $e$, proton/electron/muon/$\text{Au}^{79+}$ masses, and vacuum permittivity $\varepsilon_0$.
- **Cancellation-Free Relativistic Kinematics:** Stores proper velocity $\mathbf{u} = \gamma \mathbf{v}$ guaranteeing $|\mathbf{v}| < c$ across 12 decades of energy ($10\text{ eV}$ to $7\text{ TeV}$) and computes kinetic energy via $E_k = \frac{m_0 |\mathbf{u}|^2}{\gamma + 1}$.
- **Relativistic Boris Pusher:** Volume-preserving, time-reversible electromagnetic integrator conserving Lorentz $\gamma$ in pure magnetic fields to machine precision ($|\Delta\gamma / \gamma_0| < 10^{-14}$ over $10\,000$ steps).
- **4th-Order Yoshida Symplectic Integrator:** Exact cube-root Forest-Ruth/Yoshida splitting preserving 4D transverse betatron phase-space area ($\det M = 1$) across FODO cells.
- **Accelerator Modes:**
  1. **Synchrotron FODO Storage Ring:** Alternating focusing/defocusing quadrupoles, synchronous ramping dipoles, sextupoles, and RF accelerating cavity.
  2. **Isochronous Relativistic Cyclotron:** Radial field shaping $B_z(r) = B_0 / \sqrt{1 - (\omega_{c,0} r / c)^2}$ preserving RF phase lock into relativistic regimes.
  3. **Classical Uniform-$B_0$ Cyclotron:** Demonstrates relativistic cyclotron phase slip as $\gamma$ grows.
- **Beam Diagnostics:** Courant-Snyder Twiss parameters $(\alpha, \beta, \gamma)$ satisfying $\beta\gamma - \alpha^2 = 1$, RMS beam envelopes $(\sigma_x, \sigma_y)$, geometric & normalized emittances $(\varepsilon, \varepsilon_n)$, Larmor radius $r_L$, magnetic rigidity $B\rho$, and Liénard synchrotron radiation power $P_\gamma$.

---

## Engineering & Lifecycle Documentation
- [Software Development Life Cycle (SDLC)](docs/SDLC.md)
- [Software Requirements Specification (SRS)](docs/SRS.md)
- [Agile Development Plan & TDD Sprint Runs](docs/AGILE_TDD_PLAN.md)
- [Software Design Document (SDD)](docs/SDD.md)

---

## Quick Start

### Clone the Repository
```bash
git clone https://github.com/Mr-Kumar-Abhishek/Relativistic-Particle-Accelerator-Simulator.git
cd Relativistic-Particle-Accelerator-Simulator
```

### Run TDD Test Suite & CLI Benchmarks
```bash
cargo test
cargo run --bin particle-sim
```

### Build & Deploy WebAssembly Bundle with `wasm-builder`
```bash
# Build dist/ (compiles wasm32-unknown-unknown --release --lib + web UI + docs)
cargo run --bin wasm-builder

# Build and force-push dist/ to the gh-pages branch
cargo run --bin wasm-builder -- --deploy
```
