//! Relativistic particle representation using proper velocity u = gamma * v.

use crate::constants::{
    ELEMENTARY_CHARGE, ELECTRON_MASS, EV_TO_JOULE, GOLD_ION_CHARGE, GOLD_ION_MASS, MEV_TO_JOULE,
    MUON_MASS, PROTON_MASS, SPEED_OF_LIGHT, SPEED_OF_LIGHT_SQ, SYNCHROTRON_PREFACTOR,
};
use crate::vec3::Vec3;

/// Supported relativistic charged particle species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParticleSpecies {
    Proton,
    Electron,
    Positron,
    Muon,
    GoldIon,
    Custom { mass_kg: f64, charge_c: f64 },
}

impl ParticleSpecies {
    /// Returns the rest mass (kg) of the species.
    #[inline]
    pub fn mass_kg(self) -> f64 {
        match self {
            Self::Proton => PROTON_MASS,
            Self::Electron | Self::Positron => ELECTRON_MASS,
            Self::Muon => MUON_MASS,
            Self::GoldIon => GOLD_ION_MASS,
            Self::Custom { mass_kg, .. } => mass_kg,
        }
    }

    /// Returns the electric charge (C) of the species.
    #[inline]
    pub fn charge_c(self) -> f64 {
        match self {
            Self::Proton | Self::Positron => ELEMENTARY_CHARGE,
            Self::Electron | Self::Muon => -ELEMENTARY_CHARGE,
            Self::GoldIon => GOLD_ION_CHARGE,
            Self::Custom { charge_c, .. } => charge_c,
        }
    }

    /// Maps an integer ID (used over the WASM C-ABI) to a standard `ParticleSpecies`.
    pub fn from_id(id: u32) -> Self {
        match id {
            0 => Self::Proton,
            1 => Self::Electron,
            2 => Self::Positron,
            3 => Self::Muon,
            4 => Self::GoldIon,
            _ => Self::Proton,
        }
    }

    /// Human-readable species symbol/name.
    pub fn name(self) -> &'static str {
        match self {
            Self::Proton => "Proton (p+)",
            Self::Electron => "Electron (e-)",
            Self::Positron => "Positron (e+)",
            Self::Muon => "Muon (mu-)",
            Self::GoldIon => "Gold Ion (Au79+)",
            Self::Custom { .. } => "Custom",
        }
    }
}

/// A relativistic charged particle tracked in 3D space with proper velocity `u = gamma * v`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Particle {
    /// 3D spatial position (m).
    pub pos: Vec3,
    /// Relativistic proper velocity `u = gamma * v = p / m_0` (m/s).
    pub u: Vec3,
    /// Rest mass `m_0` (kg).
    pub mass: f64,
    /// Electric charge `q` (C).
    pub charge: f64,
    /// Active status flag (`false` when extracted or lost at aperture).
    pub alive: bool,
    /// Accumulated arc path length `s` (m).
    pub s: f64,
}

impl Particle {
    /// Constructs a particle directly from proper velocity `u = gamma * v`.
    pub fn new(species: ParticleSpecies, pos: Vec3, u: Vec3) -> Self {
        Self {
            pos,
            u,
            mass: species.mass_kg(),
            charge: species.charge_c(),
            alive: true,
            s: 0.0,
        }
    }

    /// Initializes a particle from kinetic energy in electron-volts (eV) along `direction`,
    /// avoiding floating-point cancellation at both non-relativistic and ultra-relativistic energies.
    pub fn from_kinetic_energy_ev(
        species: ParticleSpecies,
        pos: Vec3,
        direction: Vec3,
        kinetic_energy_ev: f64,
    ) -> Self {
        let mass = species.mass_kg();
        let charge = species.charge_c();
        let ek_j = kinetic_energy_ev.max(0.0) * EV_TO_JOULE;
        let rest_j = mass * SPEED_OF_LIGHT_SQ;
        let delta = ek_j / rest_j; // gamma - 1
        let u_mag = SPEED_OF_LIGHT * (delta * (delta + 2.0)).sqrt();
        let dir = direction.normalized();
        Self {
            pos,
            u: dir * u_mag,
            mass,
            charge,
            alive: true,
            s: 0.0,
        }
    }

    /// Initializes a particle from kinetic energy in MeV.
    pub fn from_kinetic_energy_mev(
        species: ParticleSpecies,
        pos: Vec3,
        direction: Vec3,
        kinetic_energy_mev: f64,
    ) -> Self {
        Self::from_kinetic_energy_ev(species, pos, direction, kinetic_energy_mev * 1.0e6)
    }

    /// Returns the Lorentz factor `gamma = sqrt(1 + |u|^2 / c^2)`.
    #[inline]
    pub fn gamma(&self) -> f64 {
        (1.0 + self.u.norm_sq() / SPEED_OF_LIGHT_SQ).sqrt()
    }

    /// Returns the dimensionless velocity vector `beta_vec = v / c = u / (gamma * c)`.
    #[inline]
    pub fn beta_vec(&self) -> Vec3 {
        self.u / (self.gamma() * SPEED_OF_LIGHT)
    }

    /// Returns the dimensionless speed `beta = v / c`.
    #[inline]
    pub fn beta(&self) -> f64 {
        self.u.norm() / (self.gamma() * SPEED_OF_LIGHT)
    }

    /// Returns the relativistic 3-velocity `v = u / gamma` (m/s), guaranteed `|v| < c`.
    #[inline]
    pub fn velocity(&self) -> Vec3 {
        self.u / self.gamma()
    }

    /// Returns the relativistic 3-momentum `p = m_0 * u = gamma * m_0 * v` (kg*m/s).
    #[inline]
    pub fn momentum(&self) -> Vec3 {
        self.u * self.mass
    }

    /// Computes kinetic energy (J) via the cancellation-free identity `E_k = m_0 * |u|^2 / (gamma + 1)`.
    #[inline]
    pub fn kinetic_energy_joules(&self) -> f64 {
        let g = self.gamma();
        self.mass * self.u.norm_sq() / (g + 1.0)
    }

    /// Returns kinetic energy in electron-volts (eV).
    #[inline]
    pub fn kinetic_energy_ev(&self) -> f64 {
        self.kinetic_energy_joules() / EV_TO_JOULE
    }

    /// Returns kinetic energy in mega-electron-volts (MeV).
    #[inline]
    pub fn kinetic_energy_mev(&self) -> f64 {
        self.kinetic_energy_joules() / MEV_TO_JOULE
    }

    /// Returns total relativistic energy `E = gamma * m_0 * c^2` (J).
    #[inline]
    pub fn total_energy_joules(&self) -> f64 {
        self.gamma() * self.mass * SPEED_OF_LIGHT_SQ
    }

    /// Returns rest mass energy `m_0 * c^2` in MeV.
    #[inline]
    pub fn rest_energy_mev(&self) -> f64 {
        (self.mass * SPEED_OF_LIGHT_SQ) / MEV_TO_JOULE
    }

    /// Returns magnetic rigidity `B*rho = p / |q| = m_0 * |u| / |q|` (T*m).
    #[inline]
    pub fn magnetic_rigidity(&self) -> f64 {
        let q = self.charge.abs();
        if q > 0.0 {
            self.mass * self.u.norm() / q
        } else {
            f64::INFINITY
        }
    }

    /// Returns the relativistic Larmor gyroradius `r_L = m_0 * |u| / (|q| * |B|)` (m) for perpendicular motion.
    #[inline]
    pub fn larmor_radius(&self, b_field_tesla: f64) -> f64 {
        let b = b_field_tesla.abs();
        let q = self.charge.abs();
        if b > 0.0 && q > 0.0 {
            (self.mass * self.u.norm()) / (q * b)
        } else {
            f64::INFINITY
        }
    }

    /// Returns the exact Larmor gyroradius using only the momentum component perpendicular to `b_vec`.
    #[inline]
    pub fn larmor_radius_in_field(&self, b_vec: Vec3) -> f64 {
        let b_norm = b_vec.norm();
        let q = self.charge.abs();
        if b_norm > 0.0 && q > 0.0 {
            let b_hat = b_vec / b_norm;
            let u_perp = self.u - b_hat * self.u.dot(b_hat);
            (self.mass * u_perp.norm()) / (q * b_norm)
        } else {
            f64::INFINITY
        }
    }

    /// Returns the relativistic angular cyclotron frequency `omega_c = |q| * |B| / (gamma * m_0)` (rad/s).
    #[inline]
    pub fn cyclotron_frequency(&self, b_field_tesla: f64) -> f64 {
        (self.charge.abs() * b_field_tesla.abs()) / (self.gamma() * self.mass)
    }

    /// Computes instantaneous Liénard synchrotron radiation power (W) for circular bending radius `rho` (m).
    #[inline]
    pub fn synchrotron_power_watts(&self, bending_radius_m: f64) -> f64 {
        if bending_radius_m <= 0.0 {
            return 0.0;
        }
        let beta_gamma_sq = self.u.norm_sq() / SPEED_OF_LIGHT_SQ;
        let beta_gamma_4 = beta_gamma_sq * beta_gamma_sq;
        (self.charge * self.charge) * SYNCHROTRON_PREFACTOR * beta_gamma_4
            / (bending_radius_m * bending_radius_m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn low_energy_avoids_catastrophic_cancellation() {
        let p = Particle::from_kinetic_energy_ev(
            ParticleSpecies::Proton,
            Vec3::zero(),
            Vec3::new(1.0, 0.0, 0.0),
            10.0,
        );
        let ek = p.kinetic_energy_ev();
        assert!((ek - 10.0).abs() / 10.0 < 1.0e-13);
        let v_classical = (2.0 * 10.0 * EV_TO_JOULE / PROTON_MASS).sqrt();
        assert!((p.velocity().norm() - v_classical).abs() / v_classical < 1.0e-7);
    }

    #[test]
    fn ultra_relativistic_lhc_proton_keeps_beta_below_one() {
        let p = Particle::from_kinetic_energy_ev(
            ParticleSpecies::Proton,
            Vec3::zero(),
            Vec3::new(0.0, 0.0, 1.0),
            7.0e12,
        );
        assert!(p.beta() < 1.0);
        assert!(p.velocity().norm() < SPEED_OF_LIGHT);
        // For kinetic energy E_k = 7 TeV, gamma = 1 + E_k / (m_p * c^2) ~= 7461.52246
        assert!((p.gamma() - 7461.522_46).abs() < 1.0e-3);
        assert!((p.kinetic_energy_ev() - 7.0e12).abs() / 7.0e12 < 1.0e-14);
    }
}
