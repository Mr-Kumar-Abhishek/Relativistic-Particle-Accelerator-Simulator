//! Fundamental physical constants (CODATA 2022 exact SI values) in double precision (`f64`).

use core::f64::consts::PI;

/// Exact speed of light in vacuum (m/s).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Squared speed of light in vacuum (m^2/s^2).
pub const SPEED_OF_LIGHT_SQ: f64 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;

/// Exact elementary charge (C).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Proton rest mass (kg), CODATA 2022.
pub const PROTON_MASS: f64 = 1.672_621_925_95e-27;

/// Electron rest mass (kg), CODATA 2022.
pub const ELECTRON_MASS: f64 = 9.109_383_713_9e-31;

/// Muon rest mass (kg), CODATA 2022.
pub const MUON_MASS: f64 = 1.883_531_627e-28;

/// Atomic mass unit (kg), CODATA 2022.
pub const ATOMIC_MASS_UNIT: f64 = 1.660_539_068_92e-27;

/// Fully stripped Gold-197 ion (Au^79+) rest mass (kg).
pub const GOLD_ION_MASS: f64 = 196.966_570_1 * ATOMIC_MASS_UNIT - 79.0 * ELECTRON_MASS;

/// Fully stripped Gold-197 ion (Au^79+) electric charge (C).
pub const GOLD_ION_CHARGE: f64 = 79.0 * ELEMENTARY_CHARGE;

/// Vacuum electric permittivity epsilon_0 (F/m), CODATA 2022.
pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_818_8e-12;

/// Electron-volt to Joule conversion factor.
pub const EV_TO_JOULE: f64 = ELEMENTARY_CHARGE;

/// Mega-electron-volt to Joule conversion factor.
pub const MEV_TO_JOULE: f64 = 1.0e6 * ELEMENTARY_CHARGE;

/// Giga-electron-volt to Joule conversion factor.
pub const GEV_TO_JOULE: f64 = 1.0e9 * ELEMENTARY_CHARGE;

/// Larmor/Lienard synchrotron radiation prefactor: q_e^2 * c / (6 * pi * eps_0) for unit charge scaling.
pub const SYNCHROTRON_PREFACTOR: f64 = SPEED_OF_LIGHT / (6.0 * PI * VACUUM_PERMITTIVITY);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codata_proton_rest_energy_mev() {
        let rest_mev = (PROTON_MASS * SPEED_OF_LIGHT_SQ) / MEV_TO_JOULE;
        assert!((rest_mev - 938.272_089_43).abs() < 1.0e-5);
    }

    #[test]
    fn codata_electron_rest_energy_mev() {
        let rest_mev = (ELECTRON_MASS * SPEED_OF_LIGHT_SQ) / MEV_TO_JOULE;
        assert!((rest_mev - 0.510_998_950_69).abs() < 1.0e-7);
    }
}
