//! Unified simulation orchestrator coordinating FODO synchrotron and cyclotron modes.

use core::f64::consts::{PI, TAU};

use crate::constants::SPEED_OF_LIGHT;
use crate::diagnostics::{BeamBunch, BeamDiagnostics};
use crate::integrator::boris_step;
use crate::lattice::{CyclotronConfig, Lattice};
use crate::particle::{Particle, ParticleSpecies};
use crate::vec3::Vec3;

/// Number of `f64` floats exported per particle in `particles_buf`.
pub const PARTICLE_STRIDE: usize = 8;
/// Number of `f64` floats exported per particle in `phase_space_buf`.
pub const PHASE_STRIDE: usize = 6;
/// Number of `f64` telemetry values exported in `diagnostics_buf`.
pub const DIAGNOSTICS_LEN: usize = 20;

/// Supported accelerator physics modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimMode {
    SynchrotronFodo = 0,
    IsochronousCyclotron = 1,
    ClassicalCyclotron = 2,
}

impl SimMode {
    pub fn from_id(id: u32) -> Self {
        match id {
            0 => Self::SynchrotronFodo,
            1 => Self::IsochronousCyclotron,
            2 => Self::ClassicalCyclotron,
            _ => Self::SynchrotronFodo,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::SynchrotronFodo => "Synchrotron FODO Storage Ring",
            Self::IsochronousCyclotron => "Isochronous Relativistic Cyclotron",
            Self::ClassicalCyclotron => "Classical Uniform-B Cyclotron",
        }
    }
}

/// High-accuracy relativistic multi-particle accelerator simulator.
pub struct AcceleratorSimulator {
    pub mode: SimMode,
    pub species: ParticleSpecies,
    pub bunch: BeamBunch,
    pub lattice: Lattice,
    pub cyclotron: CyclotronConfig,
    pub dt: f64,
    pub sim_time_s: f64,
    pub initial_rigidity_tm: f64,
    pub base_dipole_b: f64,
    pub reference_ek_mev: f64,
    pub turns_completed: f64,
    pub diagnostics: BeamDiagnostics,
    /// Flat particle buffer `[x, y, z, vx, vy, vz, gamma, ek_mev]` for zero-copy WASM access.
    pub particles_buf: Vec<f64>,
    /// Flat phase space buffer `[x_m, xp_rad, y_m, yp_rad, dphi_rad, de_mev]`.
    pub phase_space_buf: Vec<f64>,
    /// Flat telemetry array of length `DIAGNOSTICS_LEN`.
    pub diagnostics_buf: Vec<f64>,
}

impl AcceleratorSimulator {
    /// Initializes a complete accelerator simulation in the requested mode.
    pub fn new(
        mode: SimMode,
        species: ParticleSpecies,
        particle_count: usize,
        energy_mev: f64,
        b_field_t: f64,
        quad_gradient_tm: f64,
        rf_voltage_kv: f64,
    ) -> Self {
        let count = particle_count.clamp(1, 1024);
        let b_field = b_field_t.abs().max(0.05);
        let energy = energy_mev.max(0.05);

        let ref_particle = Particle::from_kinetic_energy_mev(
            species,
            Vec3::zero(),
            Vec3::new(0.0, 1.0, 0.0),
            energy,
        );
        let initial_rigidity = ref_particle.magnetic_rigidity();

        let (bunch, lattice, cyclotron, dt, ring_radius) = match mode {
            SimMode::SynchrotronFodo => {
                // Design ring radius derived directly from magnetic rigidity: R = B*rho / B0
                let r_design = (initial_rigidity / b_field).clamp(0.5, 250.0);
                let circumference = TAU * r_design;
                let num_cells = 8usize;
                let cell_len = circumference / (num_cells as f64);
                let v_rev = ref_particle.velocity().norm().max(1.0e5);
                let f_rev = v_rev / circumference;
                // Harmonic number h = 1 RF cavity synchronized to revolution frequency
                let rf_freq = f_rev;
                let lattice = Lattice::fodo_ring(
                    num_cells,
                    cell_len,
                    b_field,
                    quad_gradient_tm,
                    rf_voltage_kv * 1.0e3,
                    rf_freq,
                    0.35,
                );
                let bunch = BeamBunch::generate_ring_bunch(
                    species,
                    count,
                    r_design,
                    energy,
                    0.004 * r_design,
                    0.003,
                    0.002 * r_design,
                    0.002,
                    0.003,
                    0x4242_CAFE_1337,
                );
                let cyclo = CyclotronConfig::new(
                    b_field,
                    rf_voltage_kv * 1.0e3,
                    species.mass_kg(),
                    species.charge_c(),
                    true,
                    r_design * 1.5,
                );
                // ~240 Boris substeps per revolution for tight spatial curvature tracking
                let t_rev = circumference / v_rev;
                let dt = t_rev / 240.0;
                (bunch, lattice, cyclo, dt, r_design)
            }
            SimMode::IsochronousCyclotron | SimMode::ClassicalCyclotron => {
                let isochronous = mode == SimMode::IsochronousCyclotron;
                let r_inj = (initial_rigidity / b_field).clamp(0.03, 0.45);
                let r_max = (r_inj * 10.0).clamp(1.2, 1.85);
                let cyclo = CyclotronConfig::new(
                    b_field,
                    (rf_voltage_kv.max(25.0)) * 1.0e3,
                    species.mass_kg(),
                    species.charge_c(),
                    isochronous,
                    r_max,
                );
                let lattice = Lattice::fodo_ring(4, 1.0, b_field, 0.0, 0.0, 1.0e6, 0.0);
                let bunch = BeamBunch::generate_cyclotron_bunch(
                    species,
                    count,
                    r_inj,
                    energy,
                    0.02 * r_inj,
                    0.004,
                    0xC1C7_0001_5EED,
                );
                // ~160 steps per cyclotron period T_c = 1 / f_c
                let t_c = 1.0 / cyclo.rf_freq_hz.max(1.0);
                let dt = t_c / 160.0;
                (bunch, lattice, cyclo, dt, r_inj)
            }
        };

        let mut sim = Self {
            mode,
            species,
            bunch,
            lattice,
            cyclotron,
            dt,
            sim_time_s: 0.0,
            initial_rigidity_tm: initial_rigidity,
            base_dipole_b: b_field,
            reference_ek_mev: energy,
            turns_completed: 0.0,
            diagnostics: BeamDiagnostics::default(),
            particles_buf: vec![0.0; count * PARTICLE_STRIDE],
            phase_space_buf: vec![0.0; count * PHASE_STRIDE],
            diagnostics_buf: vec![0.0; DIAGNOSTICS_LEN],
        };
        let _ = ring_radius;
        sim.refresh_diagnostics_and_buffers();
        sim
    }

    /// Advances the entire particle bunch by `substeps` Boris steps.
    pub fn step(&mut self, substeps: usize) {
        let q_sign = self.species.charge_c().signum();
        let ap_r = match self.mode {
            SimMode::SynchrotronFodo => self.lattice.ring_radius,
            _ => self.cyclotron.max_radius_m,
        };

        for _ in 0..substeps {
            // In a Synchrotron, bending dipoles ramp synchronously with mean beam rigidity
            // so accelerated bunches stay centered on the design radius `R_0`.
            let dipole_scale = if self.mode == SimMode::SynchrotronFodo {
                let mean_rigidity = self.current_mean_rigidity();
                (mean_rigidity / self.initial_rigidity_tm.max(1.0e-12)).clamp(0.2, 100.0)
            } else {
                1.0
            };

            let t_mid = self.sim_time_s + 0.5 * self.dt;
            for p in &mut self.bunch.particles {
                if !p.alive {
                    continue;
                }
                let r_pre = (p.pos.x * p.pos.x + p.pos.y * p.pos.y).sqrt();
                if matches!(
                    self.mode,
                    SimMode::IsochronousCyclotron | SimMode::ClassicalCyclotron
                ) && r_pre >= ap_r
                {
                    // Reached extraction septum: park extracted particle at peak kinetic energy
                    continue;
                }

                let field = match self.mode {
                    SimMode::SynchrotronFodo => {
                        self.lattice.field_at_3d(p.pos, t_mid, q_sign, dipole_scale)
                    }
                    SimMode::IsochronousCyclotron | SimMode::ClassicalCyclotron => {
                        self.cyclotron.field_at(p.pos, t_mid, q_sign)
                    }
                };
                boris_step(p, field, self.dt);

                // Check physical vacuum pipe / vertical Dee plate aperture boundary
                let r = (p.pos.x * p.pos.x + p.pos.y * p.pos.y).sqrt();
                match self.mode {
                    SimMode::SynchrotronFodo => {
                        if (r - ap_r).abs() > 0.45 * ap_r || p.pos.z.abs() > 0.35 * ap_r {
                            p.alive = false;
                        }
                    }
                    SimMode::IsochronousCyclotron | SimMode::ClassicalCyclotron => {
                        if p.pos.z.abs() > ap_r * 0.4 {
                            p.alive = false;
                        }
                    }
                }
            }
            self.sim_time_s += self.dt;
        }

        self.refresh_diagnostics_and_buffers();
    }

    fn current_mean_rigidity(&self) -> f64 {
        let mut sum = 0.0;
        let mut n = 0usize;
        for p in &self.bunch.particles {
            if p.alive {
                sum += p.magnetic_rigidity();
                n += 1;
            }
        }
        if n > 0 {
            sum / (n as f64)
        } else {
            self.initial_rigidity_tm
        }
    }

    /// Recomputes statistical Twiss diagnostics and populates flat zero-copy export buffers.
    pub fn refresh_diagnostics_and_buffers(&mut self) {
        let ref_radius = match self.mode {
            SimMode::SynchrotronFodo => self.lattice.ring_radius,
            _ => {
                let mut r_sum = 0.0;
                let mut n = 0usize;
                for p in &self.bunch.particles {
                    if p.alive {
                        r_sum += (p.pos.x * p.pos.x + p.pos.y * p.pos.y).sqrt();
                        n += 1;
                    }
                }
                if n > 0 {
                    r_sum / (n as f64)
                } else {
                    0.1
                }
            }
        };

        let b_ref = match self.mode {
            SimMode::SynchrotronFodo => {
                self.base_dipole_b
                    * (self.current_mean_rigidity() / self.initial_rigidity_tm.max(1.0e-12))
            }
            _ => self.cyclotron.b0_tesla,
        };

        self.diagnostics = BeamDiagnostics::compute(&self.bunch, ref_radius, b_ref);

        let mean_ek = self.diagnostics.mean_ek_mev;
        let mut s_sum = 0.0;
        let mut alive_n = 0usize;

        for (i, p) in self.bunch.particles.iter().enumerate() {
            let v = p.velocity();
            let g = p.gamma();
            let ek_mev = p.kinetic_energy_mev();
            let p_base = i * PARTICLE_STRIDE;
            self.particles_buf[p_base] = p.pos.x;
            self.particles_buf[p_base + 1] = p.pos.y;
            self.particles_buf[p_base + 2] = p.pos.z;
            self.particles_buf[p_base + 3] = v.x;
            self.particles_buf[p_base + 4] = v.y;
            self.particles_buf[p_base + 5] = v.z;
            self.particles_buf[p_base + 6] = g;
            self.particles_buf[p_base + 7] = ek_mev;

            let r = (p.pos.x * p.pos.x + p.pos.y * p.pos.y).sqrt().max(1.0e-12);
            let cos_t = p.pos.x / r;
            let sin_t = p.pos.y / r;
            let ur = p.u.x * cos_t + p.u.y * sin_t;
            let ut = (-p.u.x * sin_t + p.u.y * cos_t).abs().max(1.0e-12);
            let theta = p.pos.y.atan2(p.pos.x);

            // Phase relative to RF frequency
            let rf_hz = match self.mode {
                SimMode::SynchrotronFodo => {
                    let v_mean = (self.diagnostics.mean_beta * SPEED_OF_LIGHT).max(1.0e4);
                    v_mean / self.lattice.circumference()
                }
                _ => self.cyclotron.rf_freq_hz,
            };
            let raw_dphi = theta - TAU * rf_hz * self.sim_time_s;
            let dphi = (raw_dphi + PI).rem_euclid(TAU) - PI;

            let ps_base = i * PHASE_STRIDE;
            self.phase_space_buf[ps_base] = r - ref_radius;
            self.phase_space_buf[ps_base + 1] = ur / ut;
            self.phase_space_buf[ps_base + 2] = p.pos.z;
            self.phase_space_buf[ps_base + 3] = p.u.z / ut;
            self.phase_space_buf[ps_base + 4] = dphi;
            self.phase_space_buf[ps_base + 5] = ek_mev - mean_ek;

            if p.alive {
                s_sum += p.s;
                alive_n += 1;
            }
        }

        let circ = match self.mode {
            SimMode::SynchrotronFodo => self.lattice.circumference(),
            _ => TAU * ref_radius.max(0.05),
        };
        self.turns_completed = if alive_n > 0 {
            (s_sum / (alive_n as f64)) / circ
        } else {
            self.turns_completed
        };

        let d = &self.diagnostics;
        self.diagnostics_buf[0] = d.mean_gamma;
        self.diagnostics_buf[1] = d.mean_beta;
        self.diagnostics_buf[2] = d.mean_ek_mev;
        self.diagnostics_buf[3] = d.mean_rigidity_tm;
        self.diagnostics_buf[4] = d.larmor_radius_m;
        self.diagnostics_buf[5] = d.cyclotron_freq_mhz;
        self.diagnostics_buf[6] = d.twiss_x.rms_pos * 1.0e3; // mm
        self.diagnostics_buf[7] = d.twiss_y.rms_pos * 1.0e3; // mm
        self.diagnostics_buf[8] = d.twiss_x.emittance_geom * 1.0e6; // mm*mrad
        self.diagnostics_buf[9] = d.twiss_x.emittance_norm * 1.0e6; // mm*mrad
        self.diagnostics_buf[10] = d.twiss_x.alpha;
        self.diagnostics_buf[11] = d.twiss_x.beta;
        self.diagnostics_buf[12] = d.twiss_x.gamma;
        self.diagnostics_buf[13] = d.energy_drift_rel;
        self.diagnostics_buf[14] = d.synchrotron_power_w;
        self.diagnostics_buf[15] = self.turns_completed;
        self.diagnostics_buf[16] = self.sim_time_s * 1.0e6; // microseconds
        self.diagnostics_buf[17] = d.survival_fraction;
        self.diagnostics_buf[18] = d.twiss_x.determinant_residual();
        self.diagnostics_buf[19] = match self.mode {
            SimMode::SynchrotronFodo => self.lattice.ring_radius,
            _ => self.cyclotron.max_radius_m,
        };
    }
}
