# Agile Development Plan & Test-Driven Development (TDD) Sprint Runs
## Project: Relativistic Particle Accelerator Simulator (`particle-accelarator`)

- **Document ID:** `AGILE-PA-001`
- **Framework:** Scrum + Strict Test-Driven Development (Red → Green → Refactor)
- **Total Planned Story Points:** `55 SP` across `5 Sprints`
- **Test Harness:** `cargo test` (`src/*` unit tests + `tests/tdd_sprints.rs` integration suite)

---

## 1. Agile Process & TDD Methodology

### 1.1 Red-Green-Refactor Discipline
Every physics feature in `particle-accelarator` is developed strictly following the TDD cycle:
1. **RED:** Write a failing unit or integration test derived directly from analytical accelerator physics equations (e.g., exact relativistic Larmor radius, Boris energy invariance, Courant-Snyder trace condition, Liouville emittance conservation).
2. **GREEN:** Implement the minimal, mathematically exact Rust code required to satisfy the physical assertion in `f64` double precision.
3. **REFACTOR:** Eliminate numerical instabilities (such as floating-point cancellation at $\beta \ll 1$), factor out shared electromagnetic helpers, verify zero compiler warnings, and re-run the entire test suite.

### 1.2 Definition of Ready (DoR) & Definition of Done (DoD)
- **Definition of Ready:** User story contains explicit physics equations, SI units, quantitative numerical error bounds (e.g., `< 1e-14`), and WASM compatibility constraints.
- **Definition of Done:**
  - Failing test written first and confirmed passing after implementation
  - `cargo test` passes 100% with zero ignored/skipped assertions
  - `cargo build --target wasm32-unknown-unknown --release --lib` succeeds
  - Logical commit pushed to GitHub (`origin`)

---

## 2. Product Backlog & Story Estimation

| Story ID | Sprint | Title & User Story | Points | Priority |
| :--- | :---: | :--- | :---: | :---: |
| **US-01** | S1 | As a physicist, I need exact CODATA 2022 constants and 3D `f64` vector math for field calculations. | 3 | P0 |
| **US-02** | S1 | As a physicist, I need proper-velocity relativistic kinematics ($\mathbf{u}=\gamma\mathbf{v}$) so $|\mathbf{v}| < c$ from $10\text{ eV}$ to $7\text{ TeV}$. | 5 | P0 |
| **US-03** | S1 | As a physicist, I need magnetic rigidity $B\rho$, Larmor radius $r_L$, cyclotron frequency $\omega_c$, and Liénard synchrotron power. | 3 | P0 |
| **US-04** | S2 | As a simulation engineer, I need a Relativistic Boris Pusher that conserves $\gamma$ to $< 10^{-14}$ over $10\,000$ turns in pure $\mathbf{B}$. | 5 | P0 |
| **US-05** | S2 | As a physicist, I need accurate cyclotron return periodicity and $\mathbf{E}\times\mathbf{B}$ drift validation. | 3 | P0 |
| **US-06** | S3 | As an optics engineer, I need `Drift`, `Dipole`, `Quadrupole`, `Sextupole`, and `RfCavity` lattice elements. | 5 | P0 |
| **US-07** | S3 | As an optics engineer, I need a 4th-order Yoshida symplectic integrator preserving transverse phase-space area. | 5 | P0 |
| **US-08** | S3 | As a physicist, I need FODO lattice stability ($f > L/4$) and instability ($f < L/4$) behavior matching Courant-Snyder theory. | 3 | P1 |
| **US-09** | S4 | As a physicist, I need RF cavity synchronous phase energy gain $\Delta E = q V_{\text{rf}} \sin\phi_s$ and longitudinal bucket oscillations. | 5 | P0 |
| **US-10** | S4 | As a physicist, I need both Classical ($B_0$) and Isochronous ($B_z = \gamma B_0$) cyclotron modes showing phase-slip elimination. | 5 | P0 |
| **US-11** | S5 | As a beam diagnostician, I need deterministic Gaussian 6D bunch generation (`XorShift64Star`) safe for WASM. | 2 | P0 |
| **US-12** | S5 | As a beam diagnostician, I need RMS emittance $(\varepsilon, \varepsilon_n)$ and Courant-Snyder Twiss $(\alpha, \beta, \gamma)$ satisfying $\beta\gamma - \alpha^2 = 1$. | 5 | P0 |
| **US-13** | S5 | As a web developer, I need zero-copy WASM FFI exports (`Float64Array` pointers) for 60 FPS browser rendering. | 3 | P0 |
| **US-14** | S5 | As a release engineer, I need a Rust `wasm-builder` CLI and interactive web dashboard deployed to `gh-pages`. | 3 | P0 |

---

## 3. Detailed TDD Sprint Run Logs

### Sprint 1: Relativistic Kinematics, CODATA Constants & 3D Vector Algebra (11 SP)
- **Sprint Goal:** Establish exact `f64` units, 3D geometry, and cancellation-free relativistic particle state representation across 12 decades of energy.
- **TDD Cycle Log:**
  1. **RED (`tests/tdd_sprints.rs::sprint1_*`):**
     - Assert $7\text{ TeV}$ LHC kinetic energy proton yields $\gamma \approx 7461.522$ and strictly $\beta < 1.0$.
     - Assert $10\text{ eV}$ non-relativistic proton matches classical $v = \sqrt{2 E_k / m_p}$ to $< 10^{-10}$ relative error and recovers exact `10.0 eV` via `kinetic_energy_ev()`.
     - Assert $B\rho = p/|q|$ and Larmor radius $r_L = p / (|q|B)$.
     - Initially fails because `constants.rs`, `vec3.rs`, `particle.rs` do not exist.
  2. **GREEN:**
     - Implemented `constants.rs` (CODATA 2022), `vec3.rs` (`Vec3` with exact vector products and rotations), and `particle.rs` (`Particle`, `ParticleSpecies`).
  3. **REFACTOR:**
     - Replaced naive `(gamma - 1.0) * m * c^2` with rationalized cancellation-free formula `m * u_sq / (gamma + 1.0)` and `u = c * ((gamma - 1.0) * (gamma + 1.0)).sqrt()`. Verified $10\text{ eV}$ accuracy improved from ~7 digits to >15 significant digits!
- **Sprint 1 Velocity:** `11 / 11 SP` (100% completed).

---

### Sprint 2: Electrodynamics & Relativistic Boris Pusher (8 SP)
- **Sprint Goal:** Implement time-reversible, volume-preserving charged particle pushing in combined 3D $\mathbf{E}$ and $\mathbf{B}$ fields.
- **TDD Cycle Log:**
  1. **RED (`tests/tdd_sprints.rs::sprint2_*`):**
     - `sprint2_boris_conserves_energy_to_machine_precision`: Inject $100\text{ MeV}$ proton into $B_z = 2.0\text{ T}$, run $10\,000$ steps ($\Delta t = 10^{-11}\text{ s}$), assert $|\gamma_{\text{final}} - \gamma_0| / \gamma_0 < 10^{-14}$.
     - `sprint2_boris_cyclotron_radius_and_period`: Assert orbital radius matches $r_L$ within $0.05\%$ and particle completes full $2\pi$ gyration after $T_c = \frac{2\pi \gamma m_0}{|q|B}$.
     - `sprint2_exb_drift_velocity`: Apply perpendicular $\mathbf{E} = (0, 10^4, 0)\text{ V/m}$ and $\mathbf{B} = (0, 0, 1.0)\text{ T}$; assert guiding-center drift matches $\mathbf{v}_d = \frac{\mathbf{E} \times \mathbf{B}}{B^2} = (10^4, 0, 0)\text{ m/s}$.
  2. **GREEN:**
     - Implemented `EmField` in `src/fields.rs` and `boris_step` in `src/integrator.rs` using half-$\mathbf{E}$ kick, $\gamma^-$ calculation, orthogonal $\mathbf{t}$ and $\mathbf{s} = 2\mathbf{t}/(1+|\mathbf{t}|^2)$ rotation, and second half-$\mathbf{E}$ kick.
  3. **REFACTOR:**
     - Extracted `Vec3::norm_sq()` inside rotation denominator to avoid unnecessary `sqrt()` calls per substep.
- **Sprint 2 Velocity:** `8 / 8 SP` (100% completed).

---

### Sprint 3: Accelerator Lattice Optics & 4th-Order Yoshida Symplectic Tracking (13 SP)
- **Sprint Goal:** Model canonical Synchrotron FODO lattices and guarantee symplectic preservation of transverse betatron phase space.
- **TDD Cycle Log:**
  1. **RED (`tests/tdd_sprints.rs::sprint3_*`):**
     - `sprint3_yoshida4_symplectic_determinant`: Compute finite-difference Jacobian matrix $M$ of 4th-order Yoshida step through quadrupole + sextupole and assert $|\det(M) - 1.0| < 10^{-11}$.
     - `sprint3_fodo_stable_vs_unstable_lattice`: Construct periodic FODO cell of total length $L = 4.0\text{ m}$. When focal length $f > L/4$, track $1\,000$ cells and assert bounded betatron amplitude ($|x| < 10\text{ mm}$). When $f < L/4$ (violating Courant-Snyder stability $|\text{Tr}(M)| \le 2$), assert exponential divergence ($|x| > 1\text{ m}$).
  2. **GREEN:**
     - Implemented `LatticeElement`, `Lattice::fodo_ring` (`src/lattice.rs`), and `yoshida4_transverse_step` (`src/integrator.rs`) using exact Yoshida weights $w_0 = -2^{1/3}/(2-2^{1/3}), w_1 = 1/(2-2^{1/3})$.
  3. **REFACTOR:**
     - Unified thin/thick multipole kick evaluations into a clean `multipole_kick` helper used across all 3 kick substeps of Yoshida-4.
- **Sprint 3 Velocity:** `13 / 13 SP` (100% completed).

---

### Sprint 4: RF Cavity Acceleration & Isochronous Cyclotron Resonance (10 SP)
- **Sprint Goal:** Enable longitudinal energy gain in synchrotron RF cavities and compare classical vs. isochronous cyclotrons at relativistic velocities.
- **TDD Cycle Log:**
  1. **RED (`tests/tdd_sprints.rs::sprint4_*`):**
     - `sprint4_rf_cavity_energy_gain`: Pass $50\text{ MeV}$ proton through `RfCavity` ($V_{\text{rf}} = 500\text{ kV}$) at $\phi_s = \pi/6$; assert energy gain equals $0.5 \times 500\text{ keV} = 250\text{ keV}$ within $1\%$.
     - `sprint4_isochronous_vs_classical_cyclotron`: Accelerate protons from $1\text{ MeV}$ in $B_0 = 1.5\text{ T}$ over $5\,000$ substeps; assert Isochronous mode ($B_z = \gamma B_0$) achieves significantly higher final kinetic energy than Classical uniform-$B_0$ mode because the classical proton slips out of RF phase resonance as $\gamma$ rises.
  2. **GREEN:**
     - Implemented `CyclotronConfig::field_at(pos, gamma, t)` and longitudinal RF phase-space tracking in `src/lattice.rs` and `src/simulator.rs`.
  3. **REFACTOR:**
     - Added phase-wrapping stabilization (`[-pi, pi]`) for longitudinal synchrotron bucket plots $(\Delta\phi, \Delta E)$.
- **Sprint 4 Velocity:** `10 / 10 SP` (100% completed).

---

### Sprint 5: Beam Diagnostics, WASM C-ABI Interface & `wasm-builder` GitHub Pages Pipeline (13 SP)
- **Sprint Goal:** Compute statistical Twiss parameters & RMS emittance, verify Liouville conservation, export zero-copy WASM FFI, and deploy live via `wasm-builder`.
- **TDD Cycle Log:**
  1. **RED (`tests/tdd_sprints.rs::sprint5_*`):**
     - `sprint5_twiss_identity_and_emittance`: Generate 128-particle Gaussian bunch, assert Courant-Snyder identity $|\beta_x \gamma_x - \alpha_x^2 - 1.0| < 10^{-12}$ and positive geometric/normalized emittance.
     - `sprint5_liouville_emittance_conservation_fodo`: Transport matched bunch across 200 linear FODO cells via Yoshida-4; assert relative emittance drift $|\varepsilon_{\text{final}} - \varepsilon_0| / \varepsilon_0 < 10^{-10}$.
     - `sprint5_wasm_ffi_buffers_valid`: Initialize simulation via `sim_init`, advance `sim_step(10)`, verify non-null pointers and valid finite telemetry metrics.
  2. **GREEN:**
     - Implemented `XorShift64Star`, `BeamBunch`, `BeamDiagnostics` (`src/diagnostics.rs`), `src/simulator.rs`, `src/wasm_api.rs`, `src/bin/wasm_builder.rs`, and `web/*`.
  3. **REFACTOR:**
     - Pre-allocated flat `Vec<f64>` output buffers inside `AcceleratorSimulator` to eliminate per-frame heap allocations in WebAssembly.
- **Sprint 5 Velocity:** `13 / 13 SP` (100% completed).

---

## 4. Sprint Burn-Up & Velocity Summary

| Sprint | Planned SP | Completed SP | Cumulative SP | Primary Verification Milestone |
| :---: | :---: | :---: | :---: | :--- |
| Sprint 1 | 11 | 11 | 11 | $10\text{ eV}$ to $7\text{ TeV}$ exact relativistic kinematics |
| Sprint 2 | 8 | 8 | 19 | Boris machine-precision $\gamma$ conservation ($<10^{-14}$) |
| Sprint 3 | 13 | 13 | 32 | Yoshida-4 symplecticity & FODO stability threshold |
| Sprint 4 | 10 | 10 | 42 | RF synchronous acceleration & Isochronous Cyclotron |
| Sprint 5 | 13 | 13 | 55 | Twiss/Liouville conservation, WASM FFI & `gh-pages` live |

**Average Velocity:** `11.0 SP / Sprint` — **100% of Product Backlog Delivered and Verified.**
