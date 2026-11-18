//! Deterministic 6D Gaussian bunch generator, Courant-Snyder Twiss parameters, and RMS emittance diagnostics.

use core::f64::consts::TAU;

use crate::constants::SPEED_OF_LIGHT;
use crate::particle::{Particle, ParticleSpecies};
use crate::vec3::Vec3;

/// Deterministic 64-bit XorShift* PRNG for reproducible WASM-compatible Gaussian bunch generation.
#[derive(Debug, Clone)]
pub struct XorShift64Star {
    state: u64,
}

impl XorShift64Star {
    /// Creates a new PRNG with a non-zero seed.
    pub fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 {
                0x9E37_79B9_7F4A_7C15
            } else {
                seed
            },
        }
    }

    /// Generates the next uniform `u64`.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Generates a uniform `f64` in `(0, 1]`.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        let v = (self.next_u64() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64));
        v.max(1.0e-15)
    }

    /// Generates a pair of independent standard normal `N(0, 1)` samples via Box-Muller transform.
    #[inline]
    pub fn next_gaussian_pair(&mut self) -> (f64, f64) {
        let u1 = self.next_f64();
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = TAU * u2;
        (r * theta.cos(), r * theta.sin())
    }
}

/// Courant-Snyder (Twiss) transverse ellipse parameters and RMS beam emittance.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TwissParameters {
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    pub rms_pos: f64,
    pub rms_angle: f64,
    pub emittance_geom: f64,
    pub emittance_norm: f64,
}

impl TwissParameters {
    /// Computes statistical Courant-Snyder Twiss parameters from paired transverse phase-space samples `(q, q')`.
    pub fn from_samples(q: &[f64], qp: &[f64], mean_beta_gamma: f64) -> Self {
        let n = q.len().min(qp.len());
        if n < 2 {
            return Self::default();
        }
        let inv_n = 1.0 / (n as f64);
        let mean_q: f64 = q[..n].iter().sum::<f64>() * inv_n;
        let mean_qp: f64 = qp[..n].iter().sum::<f64>() * inv_n;

        let mut s11 = 0.0;
        let mut s22 = 0.0;
        let mut s12 = 0.0;
        for i in 0..n {
            let dq = q[i] - mean_q;
            let dqp = qp[i] - mean_qp;
            s11 += dq * dq;
            s22 += dqp * dqp;
            s12 += dq * dqp;
        }
        s11 *= inv_n;
        s22 *= inv_n;
        s12 *= inv_n;

        let det = (s11 * s22 - s12 * s12).max(0.0);
        let emittance_geom = det.sqrt();
        if emittance_geom <= 1.0e-30 {
            return Self {
                rms_pos: s11.sqrt(),
                rms_angle: s22.sqrt(),
                ..Self::default()
            };
        }
        let beta = s11 / emittance_geom;
        let gamma = s22 / emittance_geom;
        let alpha = -s12 / emittance_geom;

        Self {
            alpha,
            beta,
            gamma,
            rms_pos: s11.sqrt(),
            rms_angle: s22.sqrt(),
            emittance_geom,
            emittance_norm: mean_beta_gamma * emittance_geom,
        }
    }

    /// Evaluates the Courant-Snyder determinant residual `|beta * gamma - alpha^2 - 1.0|`.
    #[inline]
    pub fn determinant_residual(&self) -> f64 {
        if self.emittance_geom <= 1.0e-30 {
            0.0
        } else {
            (self.beta * self.gamma - self.alpha * self.alpha - 1.0).abs()
        }
    }
}

/// Collection of macro-particles forming a relativistic charged particle bunch.
#[derive(Debug, Clone)]
pub struct BeamBunch {
    pub species: ParticleSpecies,
    pub particles: Vec<Particle>,
    pub initial_mean_gamma: f64,
}

impl BeamBunch {
    /// Generates a Gaussian 6D macro-particle bunch injected at radius `ring_radius` (m)
    /// moving counter-clockwise along `+Y`.
    pub fn generate_ring_bunch(
        species: ParticleSpecies,
        count: usize,
        ring_radius: f64,
        kinetic_energy_mev: f64,
        sigma_r_m: f64,
        sigma_rp_rad: f64,
        sigma_z_m: f64,
        sigma_zp_rad: f64,
        rel_energy_spread: f64,
        seed: u64,
    ) -> Self {
        let mut rng = XorShift64Star::new(seed);
        let mut particles = Vec::with_capacity(count);
        let ref_p = Particle::from_kinetic_energy_mev(
            species,
            Vec3::new(ring_radius, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            kinetic_energy_mev,
        );
        let u_ref = ref_p.u.norm();

        let mut gamma_sum = 0.0;
        for i in 0..count {
            let (nr, nrp) = rng.next_gaussian_pair();
            let (nz, nzp) = rng.next_gaussian_pair();
            let (ndp, nphi) = rng.next_gaussian_pair();

            // Spread slightly in azimuth around injection point for visible bunch length
            let dtheta = if i == 0 { 0.0 } else { 0.008 * nphi };
            let r = ring_radius + if i == 0 { 0.0 } else { sigma_r_m * nr };
            let z = if i == 0 { 0.0 } else { sigma_z_m * nz };
            let rp = if i == 0 { 0.0 } else { sigma_rp_rad * nrp };
            let zp = if i == 0 { 0.0 } else { sigma_zp_rad * nzp };
            let dp = if i == 0 { 0.0 } else { rel_energy_spread * ndp };

            let (sin_t, cos_t) = dtheta.sin_cos();
            let pos = Vec3::new(r * cos_t, r * sin_t, z);

            let u_mag = u_ref * (1.0 + dp).max(0.05);
            let e_r = Vec3::new(cos_t, sin_t, 0.0);
            let e_theta = Vec3::new(-sin_t, cos_t, 0.0);
            let dir = (e_theta + e_r * rp + Vec3::new(0.0, 0.0, zp)).normalized();
            let p = Particle::new(species, pos, dir * u_mag);
            gamma_sum += p.gamma();
            particles.push(p);
        }

        let initial_mean_gamma = if count > 0 {
            gamma_sum / (count as f64)
        } else {
            1.0
        };

        Self {
            species,
            particles,
            initial_mean_gamma,
        }
    }

    /// Generates a compact Gaussian bunch near origin `(r_inj, 0, 0)` for cyclotron acceleration.
    pub fn generate_cyclotron_bunch(
        species: ParticleSpecies,
        count: usize,
        r_inj_m: f64,
        kinetic_energy_mev: f64,
        sigma_r_m: f64,
        sigma_rp_rad: f64,
        seed: u64,
    ) -> Self {
        Self::generate_ring_bunch(
            species,
            count,
            r_inj_m,
            kinetic_energy_mev,
            sigma_r_m,
            sigma_rp_rad,
            0.2 * sigma_r_m,
            0.2 * sigma_rp_rad,
            0.002,
            seed,
        )
    }
}

/// Aggregated beam ensemble diagnostics snapshot.
#[derive(Debug, Clone, Copy, Default)]
pub struct BeamDiagnostics {
    pub alive_count: usize,
    pub total_count: usize,
    pub mean_gamma: f64,
    pub mean_beta: f64,
    pub mean_ek_mev: f64,
    pub mean_rigidity_tm: f64,
    pub larmor_radius_m: f64,
    pub cyclotron_freq_mhz: f64,
    pub twiss_x: TwissParameters,
    pub twiss_y: TwissParameters,
    pub energy_drift_rel: f64,
    pub synchrotron_power_w: f64,
    pub survival_fraction: f64,
}

impl BeamDiagnostics {
    /// Computes comprehensive beam diagnostics from the current bunch state.
    pub fn compute(bunch: &BeamBunch, reference_radius_m: f64, b_ref_tesla: f64) -> Self {
        let total_count = bunch.particles.len();
        let mut x_samples = Vec::with_capacity(total_count);
        let mut xp_samples = Vec::with_capacity(total_count);
        let mut y_samples = Vec::with_capacity(total_count);
        let mut yp_samples = Vec::with_capacity(total_count);

        let mut gamma_sum = 0.0;
        let mut beta_sum = 0.0;
        let mut ek_sum = 0.0;
        let mut rigidity_sum = 0.0;
        let mut larmor_sum = 0.0;
        let mut fc_sum = 0.0;
        let mut power_sum = 0.0;
        let mut bg_sum = 0.0;

        for p in &bunch.particles {
            if !p.alive {
                continue;
            }
            let r = (p.pos.x * p.pos.x + p.pos.y * p.pos.y).sqrt().max(1.0e-12);
            let cos_t = p.pos.x / r;
            let sin_t = p.pos.y / r;
            let u_mag = p.u.norm().max(1.0e-12);
            let ur = p.u.x * cos_t + p.u.y * sin_t;
            let ut = -p.u.x * sin_t + p.u.y * cos_t;

            let x_local = r - reference_radius_m;
            let xp_local = ur / ut.abs().max(1.0e-12);
            let y_local = p.pos.z;
            let yp_local = p.u.z / ut.abs().max(1.0e-12);

            x_samples.push(x_local);
            xp_samples.push(xp_local);
            y_samples.push(y_local);
            yp_samples.push(yp_local);

            let g = p.gamma();
            let b = p.beta();
            gamma_sum += g;
            beta_sum += b;
            ek_sum += p.kinetic_energy_mev();
            rigidity_sum += p.magnetic_rigidity();
            larmor_sum += p.larmor_radius(b_ref_tesla);
            fc_sum += p.cyclotron_frequency(b_ref_tesla) / TAU * 1.0e-6;
            power_sum += p.synchrotron_power_watts(reference_radius_m.max(0.05));
            bg_sum += u_mag / SPEED_OF_LIGHT;
        }

        let alive_count = x_samples.len();
        if alive_count == 0 {
            return Self {
                total_count,
                ..Self::default()
            };
        }

        let inv_n = 1.0 / (alive_count as f64);
        let mean_gamma = gamma_sum * inv_n;
        let mean_bg = bg_sum * inv_n;

        let twiss_x = TwissParameters::from_samples(&x_samples, &xp_samples, mean_bg);
        let twiss_y = TwissParameters::from_samples(&y_samples, &yp_samples, mean_bg);
        let energy_drift_rel = if bunch.initial_mean_gamma > 0.0 {
            (mean_gamma - bunch.initial_mean_gamma).abs() / bunch.initial_mean_gamma
        } else {
            0.0
        };

        Self {
            alive_count,
            total_count,
            mean_gamma,
            mean_beta: beta_sum * inv_n,
            mean_ek_mev: ek_sum * inv_n,
            mean_rigidity_tm: rigidity_sum * inv_n,
            larmor_radius_m: larmor_sum * inv_n,
            cyclotron_freq_mhz: fc_sum * inv_n,
            twiss_x,
            twiss_y,
            energy_drift_rel,
            synchrotron_power_w: power_sum * inv_n,
            survival_fraction: (alive_count as f64) / (total_count.max(1) as f64),
        }
    }
}
