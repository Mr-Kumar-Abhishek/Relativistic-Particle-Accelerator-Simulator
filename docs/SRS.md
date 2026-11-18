# Software Requirements Specification (SRS)
## Relativistic Particle Accelerator Simulator (`particle-accelarator`)

- **Document ID:** `SRS-PA-001`
- **Version:** `1.0.0`
- **Standard Alignment:** IEEE Std 830-1998 / ISO/IEC/IEEE 29148
- **Target Language & Runtime:** Rust (Edition 2021), Native CLI & WebAssembly (`wasm32-unknown-unknown`)

---

## 1. Introduction

### 1.1 Purpose
This Software Requirements Specification (SRS) defines the functional, physical, numerical accuracy, interface, and non-functional requirements for `particle-accelarator`, a high-accuracy, double-precision (`f64`) relativistic charged-particle accelerator simulator implemented in Rust and deployed to the web via WebAssembly (`wasm-builder` + GitHub Pages).

### 1.2 Scope
The system simulates single-particle and multi-particle macro-bunches inside canonical accelerator configurations:
1. **Synchrotron / FODO Storage Ring:** Bending dipole magnets, alternating-gradient focusing/defocusing quadrupoles (FODO cells), optional sextupoles, and longitudinal Radio-Frequency (RF) accelerating cavities.
2. **Isochronous & Classical Cyclotrons:** Vertical guiding magnetic fields (uniform $B_0$ or relativistic isochronous $B(r) = \gamma(r) B_0$) combined with alternating RF Dee-gap accelerating electric fields.
3. **Crossed Electromagnetic Field Validation Testbeds:** Uniform $\mathbf{B}$ gyration and $\mathbf{E} \times \mathbf{B}$ drift benchmarks.

The software delivers:
- A reusable, zero-dependency Rust physics library (`particle_accelarator`)
- A native CLI validation and benchmark runner (`particle-sim`)
- A WebAssembly FFI export layer (`wasm32-unknown-unknown`)
- A custom Rust WASM builder & GitHub Pages deployer (`wasm-builder`)
- An interactive HTML5 Canvas web dashboard displaying real-time orbits, Courant-Snyder phase space, RF buckets, conservation telemetry, and engineering docs.

### 1.3 Definitions, Acronyms, and Abbreviations
| Term | Definition |
| :--- | :--- |
| **CODATA 2022** | Committee on Data of the International Science Council recommended fundamental physical constants |
| **Proper Velocity ($\mathbf{u}$)** | Relativistic spatial four-velocity $\mathbf{u} = \gamma \mathbf{v} = \mathbf{p}/m_0$, guaranteeing $|\mathbf{v}| < c$ |
| **Lorentz Factor ($\gamma$)** | $\gamma = (1 - \beta^2)^{-1/2} = \sqrt{1 + |\mathbf{u}|^2/c^2}$ |
| **Boris Pusher** | Volume-preserving, time-reversible relativistic electromagnetic leapfrog integrator |
| **Yoshida-4** | 4th-order symplectic splitting integrator preserving Hamiltonian Poincaré invariants |
| **FODO Cell** | Periodic accelerator optics lattice: Focusing Quad – Drift – Defocusing Quad – Drift |
| **Twiss Parameters** | Courant-Snyder beam ellipse parameters $(\alpha, \beta, \gamma)$ satisfying $\beta\gamma - \alpha^2 = 1$ |
| **RMS Emittance ($\varepsilon$)** | Root-mean-square transverse phase-space area metric $\sqrt{\sigma_x^2 \sigma_{x'}^2 - \sigma_{xx'}^2}$ |
| **Magnetic Rigidity ($B\rho$)** | Beam stiffness metric $B\rho = p/|q| = m_0 |\mathbf{u}| / |q|$ in $\text{T}\cdot\text{m}$ |

---

## 2. Overall Description

### 2.1 Product Perspective
`particle-accelarator` is a standalone, self-contained simulation suite with zero external third-party runtime crate dependencies. The exact same `f64` physics core executes natively (for unit/integration tests and CLI studies) and inside browser WebAssembly linear memory.

### 2.2 User Classes and Characteristics
- **Physics Students & Researchers:** Need exact relativistic kinematics, verifiable energy conservation, and live visualization of betatron stability ($f > L/4$) vs. instability ($f < L/4$) and classical cyclotron phase slip.
- **Accelerator Engineers & Developers:** Require clean modular Rust APIs, Courant-Snyder Twiss telemetry, and reproducible TDD test suites.
- **Web Users:** Access the simulator directly on GitHub Pages with real-time parameter tuning at 60 FPS without installing Rust.

### 2.3 Operating Environment
- **Native:** Rust `rustc` / `cargo` 1.75+ on Windows, Linux, or macOS (`x86_64`).
- **WebAssembly:** `wasm32-unknown-unknown` target executing in modern ES2020+ browsers (Chrome, Firefox, Edge, Safari) with IEEE-754 `f64` support.

### 2.4 Design and Implementation Constraints
- **Zero External Dependencies:** Core library and WASM binary must rely exclusively on Rust `std` to ensure deterministic compilation, minimal `.wasm` binary size, and zero supply-chain bloat.
- **Double Precision (`f64`):** All physical state vectors, fields, and diagnostics must use 64-bit IEEE-754 floating-point arithmetic (`f64`), never `f32`.

---

## 3. Physical & Mathematical Accuracy Requirements

### 3.1 Fundamental Physical Constants (CODATA 2022)
All constants shall match SI / CODATA 2022 exact values:
- Speed of light in vacuum: $c = 299\,792\,458.0\text{ m/s}$
- Elementary charge: $e = 1.602\,176\,634 \times 10^{-19}\text{ C}$
- Proton rest mass: $m_p = 1.672\,621\,925\,95 \times 10^{-27}\text{ kg}$
- Electron rest mass: $m_e = 9.109\,383\,713\,9 \times 10^{-31}\text{ kg}$
- Muon rest mass: $m_\mu = 1.883\,531\,627 \times 10^{-28}\text{ kg}$
- Vacuum electric permittivity: $\varepsilon_0 = 8.854\,187\,818\,8 \times 10^{-12}\text{ F/m}$

### 3.2 Numerical Stability & Exact Relativistic Formulations
1. **Proper Velocity State:** Particles store $\mathbf{u} = \gamma \mathbf{v}$ instead of $\mathbf{v}$, ensuring $|\mathbf{v}| = |\mathbf{u}| / \sqrt{1 + |\mathbf{u}|^2/c^2} < c$ unconditionally, even at $100\text{ TeV}$.
2. **Cancellation-Free Kinetic Energy:** Computing $E_k = (\gamma - 1) m_0 c^2$ directly loses precision when $|\mathbf{u}| \ll c$. The simulator shall compute:
   $$E_k = \frac{m_0 |\mathbf{u}|^2}{\gamma + 1}$$
   and initialize $|\mathbf{u}|$ from kinetic energy $E_k$ via:
   $$|\mathbf{u}| = c \sqrt{(\gamma - 1)(\gamma + 1)}, \quad \text{where } \gamma - 1 = \frac{E_k}{m_0 c^2}$$
3. **Relativistic Boris Pusher:** Energy $\gamma$ in any static pure magnetic field ($\mathbf{E} = \mathbf{0}$) shall be conserved to machine precision ($|\gamma_N - \gamma_0| / \gamma_0 < 10^{-14}$ over $10\,000$ steps).
4. **4th-Order Yoshida Symplectic Integrator:** Transverse betatron tracking through magnetic multipoles shall preserve the symplectic Jacobian determinant ($\det M = 1$ within $< 10^{-12}$).

---

## 4. Functional Requirements

| Req ID | Priority | Requirement Description | Verification Method |
| :--- | :--- | :--- | :--- |
| **FR-01** | Must | Provide exact CODATA 2022 physical constants (`f64`) and $\text{eV}/\text{MeV}/\text{GeV}$ unit conversion factors. | Unit test `constants::tests` |
| **FR-02** | Must | Implement 3D Euclidean vector `Vec3` supporting dot, cross, norm, normalization, and axial rotations (`rotate_y`, `rotate_z`). | Unit test `vec3::tests` |
| **FR-03** | Must | Support multiple particle species: Proton ($p^+$), Electron ($e^-$), Positron ($e^+$), Muon ($\mu^-$), Gold Ion ($\text{Au}^{79+}$), and Custom $(m_0, q)$. | TDD Sprint 1 |
| **FR-04** | Must | Initialize particles from kinetic energy in $\text{eV}$ accurately from $10\text{ eV}$ to $7\text{ TeV}$ without cancellation or $|\mathbf{v}| \ge c$. | TDD Sprint 1 |
| **FR-05** | Must | Compute relativistic diagnostics per particle: $\gamma$, $\beta$, $\mathbf{p}$, $E_k$, $E_{\text{tot}}$, magnetic rigidity $B\rho$, Larmor radius $r_L$, cyclotron frequency $\omega_c$, and Liénard synchrotron power $P_\gamma$. | TDD Sprint 1 |
| **FR-06** | Must | Implement the Relativistic Boris Pusher for arbitrary 3D $\mathbf{E}(\mathbf{r}, t), \mathbf{B}(\mathbf{r}, t)$ fields. | TDD Sprint 2 |
| **FR-07** | Must | Guarantee exact energy conservation ($|\Delta\gamma/\gamma_0| < 10^{-14}$) over $10\,000+$ cyclotron gyrations and accurate $\mathbf{E} \times \mathbf{B}$ drift velocity. | TDD Sprint 2 |
| **FR-08** | Must | Model standard accelerator lattice elements: `Drift`, `Dipole`, `Quadrupole`, `Sextupole`, and `RfCavity`. | TDD Sprint 3 |
| **FR-09** | Must | Construct periodic FODO storage ring lattices and evaluate local electromagnetic fields along orbital arc length $s$. | TDD Sprint 3 |
| **FR-10** | Must | Implement a 4th-order Yoshida symplectic integrator for transverse betatron phase space $(x, p_x, y, p_y)$ preserving phase-space area. | TDD Sprint 3 |
| **FR-11** | Must | Demonstrate Courant-Snyder betatron stability ($f > L/4$) and instability ($f < L/4$) in FODO lattices. | TDD Sprint 3 |
| **FR-12** | Must | Model longitudinal RF cavity acceleration $\Delta E = q V_{\text{rf}} \sin(\phi_s)$ and synchrotron phase oscillations $(\Delta\phi, \Delta E)$. | TDD Sprint 4 |
| **FR-13** | Must | Simulate both Classical uniform-$B_0$ Cyclotrons and Isochronous Cyclotrons ($B_z(r) = \gamma(r) B_0$), demonstrating phase-lock preservation at relativistic energies. | TDD Sprint 4 |
| **FR-14** | Must | Generate reproducible 6D Gaussian macro-particle bunches using a deterministic `XorShift64Star` Box-Muller PRNG (no OS RNG syscalls). | TDD Sprint 5 |
| **FR-15** | Must | Compute bunch Courant-Snyder Twiss parameters $(\alpha, \beta, \gamma)$, RMS beam envelopes $(\sigma_x, \sigma_y)$, geometric emittance $\varepsilon$, and normalized emittance $\varepsilon_n$, satisfying $\beta\gamma - \alpha^2 = 1$. | TDD Sprint 5 |
| **FR-16** | Must | Export a zero-copy C-ABI WebAssembly interface (`sim_init`, `sim_step`, `sim_get_particles_ptr`, `sim_get_phase_space_ptr`, `sim_get_diagnostics_ptr`, `sim_get_particle_count`). | TDD Sprint 5 |
| **FR-17** | Must | Provide a Rust `wasm-builder` CLI tool that compiles `wasm32-unknown-unknown` in release mode, validates `.wasm` magic bytes, bundles `web/` + `docs/` into `dist/`, and deploys to `gh-pages`. | Step 4 & 5 Verification |
| **FR-18** | Must | Provide an interactive browser UI visualizing orbits, transverse Twiss ellipses, longitudinal RF phase space, conservation plots, and embedded SDLC/SRS/Agile/SDD documentation. | Step 4 & 5 Verification |

---

## 5. External Interface Requirements

### 5.1 Native CLI Interface (`particle-sim`)
- Outputs formatted physics verification summaries for:
  - Synchrotron FODO Proton Ring ($450\text{ GeV}$ / adjustable)
  - Isochronous vs. Classical Cyclotron Acceleration Comparison
  - 10,000-step Boris Energy Conservation Benchmark

### 5.2 WebAssembly C-ABI Interface (`particle_accelarator.wasm`)
- `sim_init(mode: u32, particle_count: u32, energy_mev: f64, b_field_t: f64, quad_gradient: f64, rf_voltage_kv: f64, isochronous: u32)`
- `sim_step(substeps: u32)`
- `sim_get_particles_ptr() -> *const f64` (`8 * N` floats: `[x, y, z, vx, vy, vz, gamma, ek_mev]`)
- `sim_get_phase_space_ptr() -> *const f64` (`6 * N` floats: `[x_m, xp_rad, y_m, yp_rad, dphi_rad, de_mev]`)
- `sim_get_diagnostics_ptr() -> *const f64` (`20` floats of real-time telemetry)
- `sim_get_particle_count() -> u32`
- `sim_get_ring_radius() -> f64`

---

## 6. Non-Functional Requirements

| Req ID | Category | Specification |
| :--- | :--- | :--- |
| **NFR-01** | Numerical Accuracy | All particle invariants in pure static $\mathbf{B}$ fields conserved to relative error $< 10^{-14}$ over $10\,000$ Boris steps. |
| **NFR-02** | Symplecticity | Courant-Snyder Twiss determinant error $|\beta_x \gamma_x - \alpha_x^2 - 1| < 10^{-12}$; linear FODO emittance drift $< 10^{-10}$. |
| **NFR-03** | Performance | Single frame step (256 particles $\times$ 8 substeps) completes in $< 4\text{ ms}$ in browser WebAssembly (well within $16.6\text{ ms}$ 60 FPS budget). |
| **NFR-04** | Determinism | Identical simulation parameters produce bit-for-bit reproducible trajectories across native `x86_64` and `wasm32` targets. |
| **NFR-05** | Footprint | Compiled release `.wasm` binary size $< 250\text{ KB}$ with zero external npm or crate runtime dependencies. |
| **NFR-06** | Portability | Builds cleanly on Windows, Linux, macOS, and GitHub Actions Ubuntu runners. |
| **NFR-07** | Maintainability | Modular single-responsibility Rust modules with docstrings and strict TDD regression suite. |
| **NFR-08** | Accessibility | Deployed publicly via GitHub Pages (`gh-pages`) with `.nojekyll` and built-in documentation viewer. |

---

## 7. Requirements Traceability Matrix (RTM)

| Requirement IDs | Design Module (`src/`) | Verification Test / Script |
| :--- | :--- | :--- |
| `FR-01`, `FR-02` | `src/constants.rs`, `src/vec3.rs` | `tests/tdd_sprints.rs::sprint1_*`, module unit tests |
| `FR-03`, `FR-04`, `FR-05` | `src/particle.rs` | `tests/tdd_sprints.rs::sprint1_*` |
| `FR-06`, `FR-07` | `src/fields.rs`, `src/integrator.rs` | `tests/tdd_sprints.rs::sprint2_*` |
| `FR-08`, `FR-09`, `FR-10`, `FR-11` | `src/lattice.rs`, `src/integrator.rs` | `tests/tdd_sprints.rs::sprint3_*` |
| `FR-12`, `FR-13` | `src/lattice.rs`, `src/simulator.rs` | `tests/tdd_sprints.rs::sprint4_*` |
| `FR-14`, `FR-15`, `FR-16` | `src/diagnostics.rs`, `src/wasm_api.rs` | `tests/tdd_sprints.rs::sprint5_*` |
| `FR-17`, `FR-18` | `src/bin/wasm_builder.rs`, `web/*` | `cargo run --bin wasm-builder`, GitHub Pages deploy |
