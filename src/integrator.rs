//! Symplectic and phase-space volume-preserving integrators for relativistic tracking.

use crate::constants::SPEED_OF_LIGHT_SQ;
use crate::fields::EmField;
use crate::particle::Particle;

/// Exact 4th-order Forest-Ruth / Yoshida splitting coefficients:
/// w0 = -2^(1/3) / (2 - 2^(1/3)), w1 = 1 / (2 - 2^(1/3)).
const CBRT_2: f64 = 1.259_921_049_894_873_2; // 2^(1/3) in full f64 precision
const W0: f64 = -CBRT_2 / (2.0 - CBRT_2);
const W1: f64 = 1.0 / (2.0 - CBRT_2);

const YOSHIDA_C: [f64; 4] = [0.5 * W1, 0.5 * (W0 + W1), 0.5 * (W0 + W1), 0.5 * W1];
const YOSHIDA_D: [f64; 3] = [W1, W0, W1];

/// Advances a relativistic particle over time step `dt` (s) in electromagnetic field `(E, B)`
/// using the volume-preserving, time-reversible Relativistic Boris Pusher.
#[inline]
pub fn boris_step(particle: &mut Particle, field: EmField, dt: f64) {
    if !particle.alive {
        return;
    }

    let q_dt_half_m = (particle.charge * dt) / (2.0 * particle.mass);

    // 1. Half electric impulse
    let u_minus = particle.u + field.e * q_dt_half_m;

    // 2. Intermediate relativistic Lorentz factor
    let gamma_minus = (1.0 + u_minus.norm_sq() / SPEED_OF_LIGHT_SQ).sqrt();

    // 3. Magnetic rotation via exact Cayley transform (orthogonal -> preserves |u| exactly)
    let t = field.b * (q_dt_half_m / gamma_minus);
    let s = t * (2.0 / (1.0 + t.norm_sq()));
    let u_prime = u_minus + u_minus.cross(t);
    let u_plus = u_minus + u_prime.cross(s);

    // 4. Second half electric impulse
    let u_next = u_plus + field.e * q_dt_half_m;
    let gamma_next = (1.0 + u_next.norm_sq() / SPEED_OF_LIGHT_SQ).sqrt();

    // 5. Update proper velocity, relativistic position drift, and path length
    let v_next = u_next / gamma_next;
    let dr = v_next * dt;
    particle.u = u_next;
    particle.pos += dr;
    particle.s += dr.norm();
}

/// Evaluates the exact curl-free magnetic multipole kick for transverse betatron coordinates.
#[inline]
fn multipole_kick(
    x: f64,
    px: &mut f64,
    y: f64,
    py: &mut f64,
    step_s: f64,
    k1: f64,
    k2: f64,
) {
    let fx = -k1 * x - 0.5 * k2 * (x * x - y * y);
    let fy = k1 * y + k2 * x * y;
    *px += fx * step_s;
    *py += fy * step_s;
}

/// Advances 4D transverse betatron phase space `(x, px, y, py)` along longitudinal step `ds` (m)
/// through a combined quadrupole (`k1`, m^-2) and sextupole (`k2`, m^-3) element using a
/// 4th-order Yoshida symplectic integrator.
pub fn yoshida4_transverse_step(
    x: &mut f64,
    px: &mut f64,
    y: &mut f64,
    py: &mut f64,
    ds: f64,
    k1: f64,
    k2: f64,
) {
    *x += *px * (YOSHIDA_C[0] * ds);
    *y += *py * (YOSHIDA_C[0] * ds);
    multipole_kick(*x, px, *y, py, YOSHIDA_D[0] * ds, k1, k2);

    *x += *px * (YOSHIDA_C[1] * ds);
    *y += *py * (YOSHIDA_C[1] * ds);
    multipole_kick(*x, px, *y, py, YOSHIDA_D[1] * ds, k1, k2);

    *x += *px * (YOSHIDA_C[2] * ds);
    *y += *py * (YOSHIDA_C[2] * ds);
    multipole_kick(*x, px, *y, py, YOSHIDA_D[2] * ds, k1, k2);

    *x += *px * (YOSHIDA_C[3] * ds);
    *y += *py * (YOSHIDA_C[3] * ds);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::ParticleSpecies;
    use crate::vec3::Vec3;

    #[test]
    fn boris_conserves_gamma_in_uniform_magnetic_field() {
        let mut p = Particle::from_kinetic_energy_mev(
            ParticleSpecies::Proton,
            Vec3::zero(),
            Vec3::new(1.0, 0.0, 0.0),
            100.0,
        );
        let gamma0 = p.gamma();
        let field = EmField::magnetic(Vec3::new(0.0, 0.0, 2.0));
        for _ in 0..10_000 {
            boris_step(&mut p, field, 1.0e-11);
        }
        let rel_err = (p.gamma() - gamma0).abs() / gamma0;
        assert!(rel_err < 1.0e-14, "relative gamma drift: {rel_err}");
    }

    #[test]
    fn yoshida_coefficients_sum_to_unity() {
        let c_sum: f64 = YOSHIDA_C.iter().sum();
        let d_sum: f64 = YOSHIDA_D.iter().sum();
        assert!((c_sum - 1.0).abs() < 1.0e-15);
        assert!((d_sum - 1.0).abs() < 1.0e-15);
    }
}
