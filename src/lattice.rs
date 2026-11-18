//! Accelerator lattice elements, periodic FODO rings, and cyclotron field models.

use core::f64::consts::{PI, TAU};

use crate::constants::SPEED_OF_LIGHT;
use crate::fields::EmField;
use crate::integrator::yoshida4_transverse_step;
use crate::vec3::Vec3;

/// Individual electromagnetic element inside an accelerator ring or beamline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LatticeElement {
    /// Field-free drift space of `length` (m).
    Drift { length: f64 },
    /// Sector bending dipole magnet of arc `length` (m), vertical field `b_y` (T), and design radius `bending_radius` (m).
    Dipole {
        length: f64,
        b_y: f64,
        bending_radius: f64,
    },
    /// Magnetic quadrupole of `length` (m) and field gradient `gradient` (T/m).
    /// Positive `gradient` focuses horizontally (`x`) and defocuses vertically (`y`).
    Quadrupole { length: f64, gradient: f64 },
    /// Magnetic sextupole of `length` (m) and strength `strength` (T/m^2).
    Sextupole { length: f64, strength: f64 },
    /// Longitudinal Radio-Frequency accelerating cavity.
    RfCavity {
        length: f64,
        voltage_v: f64,
        frequency_hz: f64,
        phase_rad: f64,
    },
}

impl LatticeElement {
    /// Returns the physical arc length (m) of this lattice element.
    #[inline]
    pub fn length(&self) -> f64 {
        match *self {
            Self::Drift { length }
            | Self::Dipole { length, .. }
            | Self::Quadrupole { length, .. }
            | Self::Sextupole { length, .. }
            | Self::RfCavity { length, .. } => length,
        }
    }
}

/// A periodic sequence of accelerator lattice elements forming a storage ring or synchrotron.
#[derive(Debug, Clone)]
pub struct Lattice {
    pub elements: Vec<LatticeElement>,
    pub ring_radius: f64,
    pub guiding_b_y: f64,
    pub dipole_fill_fraction: f64,
    circumference_cached: f64,
}

impl Lattice {
    /// Constructs a lattice from an ordered list of elements and design ring radius (m).
    pub fn new(elements: Vec<LatticeElement>, ring_radius: f64) -> Self {
        let mut total_len = 0.0;
        let mut dipole_len = 0.0;
        let mut dipole_b_sum = 0.0;
        for elem in &elements {
            let l = elem.length();
            total_len += l;
            if let LatticeElement::Dipole { b_y, .. } = *elem {
                dipole_len += l;
                dipole_b_sum += b_y * l;
            }
        }
        let dipole_fill_fraction = if total_len > 0.0 && dipole_len > 0.0 {
            dipole_len / total_len
        } else {
            1.0
        };
        let guiding_b_y = if total_len > 0.0 && dipole_len > 0.0 {
            dipole_b_sum / total_len
        } else {
            0.0
        };
        Self {
            elements,
            ring_radius: ring_radius.max(1.0e-3),
            guiding_b_y,
            dipole_fill_fraction,
            circumference_cached: total_len.max(1.0e-6),
        }
    }

    /// Constructs a periodic FODO storage ring with `num_cells` cells, each of length `cell_length` (m).
    pub fn fodo_ring(
        num_cells: usize,
        cell_length: f64,
        dipole_b_y: f64,
        quad_gradient: f64,
        rf_voltage_v: f64,
        rf_freq_hz: f64,
        rf_phase_rad: f64,
    ) -> Self {
        let cells = num_cells.max(1);
        let total_circumference = (cells as f64) * cell_length;
        let ring_radius = total_circumference / TAU;

        let l_quad = 0.10 * cell_length;
        let l_dipole = 0.30 * cell_length;
        let l_gap = 0.10 * cell_length;

        let mut elements = Vec::with_capacity(cells * 6);
        for cell_idx in 0..cells {
            elements.push(LatticeElement::Quadrupole {
                length: l_quad,
                gradient: quad_gradient,
            });
            elements.push(LatticeElement::Dipole {
                length: l_dipole,
                b_y: dipole_b_y,
                bending_radius: ring_radius,
            });
            if cell_idx == 0 && rf_voltage_v.abs() > 0.0 {
                elements.push(LatticeElement::RfCavity {
                    length: l_gap,
                    voltage_v: rf_voltage_v,
                    frequency_hz: rf_freq_hz,
                    phase_rad: rf_phase_rad,
                });
            } else {
                elements.push(LatticeElement::Drift { length: l_gap });
            }
            elements.push(LatticeElement::Quadrupole {
                length: l_quad,
                gradient: -quad_gradient,
            });
            elements.push(LatticeElement::Dipole {
                length: l_dipole,
                b_y: dipole_b_y,
                bending_radius: ring_radius,
            });
            elements.push(LatticeElement::Drift { length: l_gap });
        }

        let mut lat = Self::new(elements, ring_radius);
        lat.guiding_b_y = dipole_b_y;
        lat
    }

    /// Total circumference (m) of the lattice.
    #[inline]
    pub fn circumference(&self) -> f64 {
        self.circumference_cached
    }

    /// Looks up the active lattice element and offset within it at longitudinal arc coordinate `s` (m).
    pub fn element_at_s(&self, s: f64) -> (usize, &LatticeElement, f64) {
        let circ = self.circumference_cached;
        let s_mod = s.rem_euclid(circ);
        let mut accum = 0.0;
        for (idx, elem) in self.elements.iter().enumerate() {
            let l = elem.length();
            if s_mod < accum + l || idx + 1 == self.elements.len() {
                return (idx, elem, s_mod - accum);
            }
            accum += l;
        }
        (0, &self.elements[0], 0.0)
    }

    /// Evaluates the local 3D electromagnetic field in the horizontal X-Y ring plane (Z vertical)
    /// for a particle at position `pos`, with particle charge polarity `charge_sign` (+1 or -1)
    /// and synchronous magnetic field ramp factor `dipole_scale`.
    pub fn field_at_3d(
        &self,
        pos: Vec3,
        t: f64,
        charge_sign: f64,
        dipole_scale: f64,
    ) -> EmField {
        let r = (pos.x * pos.x + pos.y * pos.y).sqrt();
        if r < 1.0e-9 {
            return EmField::zero();
        }
        let theta = pos.y.atan2(pos.x).rem_euclid(TAU);
        let s = theta * self.ring_radius;
        let (_, elem, _) = self.element_at_s(s);

        let e_r = Vec3::new(pos.x / r, pos.y / r, 0.0);
        let e_theta = Vec3::new(-pos.y / r, pos.x / r, 0.0);
        let dr = r - self.ring_radius;
        let z_vert = pos.z;

        // For counter-clockwise motion (e_theta) of positive charge, inward centripetal force requires -Z B-field.
        let pol = if charge_sign >= 0.0 { 1.0 } else { -1.0 };
        let b_guide = self.guiding_b_y * dipole_scale;
        let bz_base = -pol * b_guide * (1.0 - 0.05 * dr / self.ring_radius);
        let br_base = -pol * b_guide * (-0.05 * z_vert / self.ring_radius);

        match *elem {
            LatticeElement::Drift { .. } | LatticeElement::Dipole { .. } => {
                EmField::magnetic(e_r * br_base + Vec3::new(0.0, 0.0, bz_base))
            }
            LatticeElement::Quadrupole { gradient, .. } => {
                let bz_quad = -pol * gradient * dr;
                let br_quad = -pol * gradient * z_vert;
                EmField::magnetic(e_r * (br_base + br_quad) + Vec3::new(0.0, 0.0, bz_base + bz_quad))
            }
            LatticeElement::Sextupole { strength, .. } => {
                let bz_sex = -pol * 0.5 * strength * (dr * dr - z_vert * z_vert);
                let br_sex = -pol * strength * dr * z_vert;
                EmField::magnetic(e_r * (br_base + br_sex) + Vec3::new(0.0, 0.0, bz_base + bz_sex))
            }
            LatticeElement::RfCavity {
                length,
                voltage_v,
                frequency_hz,
                phase_rad,
            } => {
                let e_mag = (voltage_v / length.max(1.0e-3))
                    * (TAU * frequency_hz * t + phase_rad).sin();
                EmField::new(
                    e_theta * (pol * e_mag),
                    e_r * br_base + Vec3::new(0.0, 0.0, bz_base),
                )
            }
        }
    }

    /// Tracks 4D transverse betatron coordinates `(x, px, y, py)` through one full pass
    /// of this lattice using the 4th-order Yoshida symplectic integrator.
    pub fn track_yoshida4_pass(
        &self,
        x: &mut f64,
        px: &mut f64,
        y: &mut f64,
        py: &mut f64,
        rigidity_tm: f64,
        substeps_per_element: usize,
    ) {
        let steps = substeps_per_element.max(1);
        for elem in &self.elements {
            let l = elem.length();
            let ds = l / (steps as f64);
            let (k1, k2) = match *elem {
                LatticeElement::Quadrupole { gradient, .. } => (gradient / rigidity_tm, 0.0),
                LatticeElement::Sextupole { strength, .. } => (0.0, strength / rigidity_tm),
                _ => (0.0, 0.0),
            };
            for _ in 0..steps {
                yoshida4_transverse_step(x, px, y, py, ds, k1, k2);
            }
        }
    }
}

/// Configuration for a Classical or Relativistic Isochronous Cyclotron.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CyclotronConfig {
    /// Central vertical magnetic flux density `B_0` (T).
    pub b0_tesla: f64,
    /// Width of the RF accelerating Dee gap along `|x| < dee_gap_m / 2` (m).
    pub dee_gap_m: f64,
    /// Peak RF gap voltage `V_rf` (V).
    pub rf_voltage_v: f64,
    /// RF drive frequency (Hz), typically matched to `|q| * B_0 / (2 * pi * m_0)`.
    pub rf_freq_hz: f64,
    /// If true, uses radial relativistic field shaping `B_z(r) = B_0 * gamma_iso(r)`.
    pub isochronous: bool,
    /// Extraction aperture radius (m).
    pub max_radius_m: f64,
}

impl CyclotronConfig {
    /// Creates a cyclotron configuration tuned to the non-relativistic cyclotron frequency of `(mass, charge)`.
    pub fn new(
        b0_tesla: f64,
        rf_voltage_v: f64,
        mass_kg: f64,
        charge_c: f64,
        isochronous: bool,
        max_radius_m: f64,
    ) -> Self {
        let rf_freq_hz = (charge_c.abs() * b0_tesla.abs()) / (TAU * mass_kg);
        Self {
            b0_tesla,
            dee_gap_m: (0.06 * max_radius_m).clamp(0.02, 0.25),
            rf_voltage_v,
            rf_freq_hz,
            isochronous,
            max_radius_m,
        }
    }

    /// Evaluates the 3D electromagnetic field `(E, B)` at spatial coordinate `pos` and time `t`.
    /// Note: For isochronous mode, `B_z(r) = B_0 / sqrt(1 - (omega_0 * r / c)^2)` depends purely on
    /// radius `r` so it remains a true static magnetostatic field with fixed orbital frequency.
    pub fn field_at(&self, pos: Vec3, t: f64, charge_sign: f64) -> EmField {
        let pol = if charge_sign >= 0.0 { 1.0 } else { -1.0 };
        let r = (pos.x * pos.x + pos.y * pos.y).sqrt();

        let bz_mag = if self.isochronous {
            let omega_0 = TAU * self.rf_freq_hz;
            let beta_r = (omega_0 * r / SPEED_OF_LIGHT).min(0.992);
            let gamma_iso = 1.0 / (1.0 - beta_r * beta_r).sqrt();
            self.b0_tesla * gamma_iso
        } else {
            self.b0_tesla
        };

        // Positive particle moving CCW: e_theta x (-e_z) = -e_r (inward centripetal force)
        let b_vec = Vec3::new(0.0, 0.0, -pol * bz_mag);

        // Alternating RF electric field inside the Dee gap (|x| < dee_gap_m / 2).
        // Particle starts at (+r, 0) moving CCW (+Y), so it crosses x = 0 moving -X at t = T_c/4
        // and moving +X at t = 3*T_c/4. -sin(omega_rf * t) peaks at -1 and +1 at those crossings.
        let half_gap = 0.5 * self.dee_gap_m;
        let e_vec = if pos.x.abs() <= half_gap {
            let taper = 0.5 * (1.0 + (PI * pos.x / half_gap).cos());
            let ex = -pol * (self.rf_voltage_v / self.dee_gap_m) * (TAU * self.rf_freq_hz * t).sin()
                * (2.0 * taper);
            Vec3::new(ex, 0.0, 0.0)
        } else {
            Vec3::zero()
        };

        EmField::new(e_vec, b_vec)
    }
}
