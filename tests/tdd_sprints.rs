//! End-to-end TDD Sprint Integration Test Suite (Sprints 1 through 5).

use core::f64::consts::PI;

use particle_accelarator::constants::{
    ELECTRON_MASS, EV_TO_JOULE, PROTON_MASS, SPEED_OF_LIGHT, SPEED_OF_LIGHT_SQ, VACUUM_PERMITTIVITY,
};
use particle_accelarator::diagnostics::{BeamBunch, BeamDiagnostics, TwissParameters};
use particle_accelarator::fields::EmField;
use particle_accelarator::integrator::{boris_step, yoshida4_transverse_step};
use particle_accelarator::lattice::{Lattice, LatticeElement};
use particle_accelarator::particle::{Particle, ParticleSpecies};
use particle_accelarator::simulator::{AcceleratorSimulator, SimMode};
use particle_accelarator::vec3::Vec3;
use particle_accelarator::wasm_api::{
    sim_get_diagnostics_ptr, sim_get_particle_count, sim_get_particles_ptr,
    sim_get_phase_space_ptr, sim_init, sim_step,
};

fn det_4x4(m: [[f64; 4]; 4]) -> f64 {
    let mut a = m;
    let mut det = 1.0;
    for col in 0..4 {
        let mut pivot = col;
        for row in (col + 1)..4 {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        if a[pivot][col].abs() < 1.0e-30 {
            return 0.0;
        }
        if pivot != col {
            a.swap(col, pivot);
            det = -det;
        }
        det *= a[col][col];
        let inv = 1.0 / a[col][col];
        for row in (col + 1)..4 {
            let factor = a[row][col] * inv;
            for k in col..4 {
                a[row][k] -= factor * a[col][k];
            }
        }
    }
    det
}

#[test]
fn sprint1_relativistic_kinematics_from_ev_to_tev() {
    // Low-energy non-relativistic limit: 10 eV proton without cancellation
    let low = Particle::from_kinetic_energy_ev(
        ParticleSpecies::Proton,
        Vec3::zero(),
        Vec3::new(1.0, 0.0, 0.0),
        10.0,
    );
    assert!((low.kinetic_energy_ev() - 10.0).abs() / 10.0 < 1.0e-13);
    let v_classical = (2.0 * 10.0 * EV_TO_JOULE / PROTON_MASS).sqrt();
    assert!((low.velocity().norm() - v_classical).abs() / v_classical < 1.0e-7);

    // Ultra-relativistic limit: 7 TeV LHC proton
    let high = Particle::from_kinetic_energy_ev(
        ParticleSpecies::Proton,
        Vec3::zero(),
        Vec3::new(0.0, 1.0, 0.0),
        7.0e12,
    );
    let expected_gamma = 1.0 + (7.0e12 * EV_TO_JOULE) / (PROTON_MASS * SPEED_OF_LIGHT_SQ);
    assert!(high.beta() < 1.0);
    assert!(high.velocity().norm() < SPEED_OF_LIGHT);
    assert!((high.gamma() - expected_gamma).abs() / expected_gamma < 1.0e-14);
}

#[test]
fn sprint1_magnetic_rigidity_larmor_and_synchrotron_power() {
    let p = Particle::from_kinetic_energy_mev(
        ParticleSpecies::Proton,
        Vec3::zero(),
        Vec3::new(1.0, 0.0, 0.0),
        500.0,
    );
    let expected_rigidity = p.momentum().norm() / p.charge.abs();
    assert!((p.magnetic_rigidity() - expected_rigidity).abs() < 1.0e-15);
    assert!((p.larmor_radius(2.0) - expected_rigidity / 2.0).abs() < 1.0e-15);

    let e = Particle::from_kinetic_energy_mev(
        ParticleSpecies::Electron,
        Vec3::zero(),
        Vec3::new(1.0, 0.0, 0.0),
        2_000.0,
    );
    let rho = 10.0;
    let bg = e.beta() * e.gamma();
    let expected_power = (EV_TO_JOULE * EV_TO_JOULE * SPEED_OF_LIGHT / (6.0 * PI * VACUUM_PERMITTIVITY))
        * bg.powi(4)
        / (rho * rho);
    let rel_err = (e.synchrotron_power_watts(rho) - expected_power).abs() / expected_power;
    assert!(rel_err < 1.0e-13, "synchrotron power rel_err: {rel_err}");
    assert!(ELECTRON_MASS < PROTON_MASS);
}

#[test]
fn sprint2_boris_conserves_energy_to_machine_precision() {
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
    let rel_drift = (p.gamma() - gamma0).abs() / gamma0;
    assert!(rel_drift < 1.0e-14, "gamma drift = {rel_drift}");
}

#[test]
fn sprint2_boris_cyclotron_radius_and_period() {
    let b_z = 1.5;
    let mut p = Particle::from_kinetic_energy_mev(
        ParticleSpecies::Proton,
        Vec3::zero(),
        Vec3::new(1.0, 0.0, 0.0),
        20.0,
    );
    let r_l = p.larmor_radius(b_z);
    let t_c = 2.0 * PI / p.cyclotron_frequency(b_z);
    let steps = 1_000;
    let dt = t_c / (steps as f64);
    let field = EmField::magnetic(Vec3::new(0.0, 0.0, b_z));

    for i in 1..=steps {
        boris_step(&mut p, field, dt);
        if i == steps / 2 {
            // Half-period diameter across gyration center (0, -r_l, 0)
            let diameter_err = (p.pos.norm() - 2.0 * r_l).abs() / (2.0 * r_l);
            assert!(diameter_err < 1.0e-4, "diameter_err = {diameter_err}");
        }
    }
    let return_err = p.pos.norm() / r_l;
    assert!(return_err < 1.0e-3, "return_err = {return_err}");
}

#[test]
fn sprint2_exb_drift_velocity() {
    let mut p = Particle::new(ParticleSpecies::Proton, Vec3::zero(), Vec3::zero());
    let e_field = Vec3::new(0.0, 1.0e4, 0.0);
    let b_field = Vec3::new(0.0, 0.0, 1.0);
    let field = EmField::new(e_field, b_field);

    let t_c = 2.0 * PI / p.cyclotron_frequency(1.0);
    let n_periods = 10;
    let steps_per_period = 400;
    let dt = t_c / (steps_per_period as f64);

    for _ in 0..(n_periods * steps_per_period) {
        boris_step(&mut p, field, dt);
    }

    let total_t = (n_periods as f64) * t_c;
    let avg_drift = p.pos / total_t;
    let expected_vd = e_field.cross(b_field) / b_field.norm_sq(); // (1e4, 0, 0)
    let err = (avg_drift - expected_vd).norm() / expected_vd.norm();
    assert!(err < 1.0e-3, "ExB drift error = {err}, avg = {:?}", avg_drift);
}

#[test]
fn sprint3_yoshida4_symplectic_determinant() {
    let base = [0.002_f64, 0.0005, -0.001, 0.0003];
    let eps = 1.0e-7;
    let ds = 0.25;
    let k1 = 1.8;
    let k2 = 12.0;

    let mut jac = [[0.0_f64; 4]; 4];
    for j in 0..4 {
        let [mut px0, mut px1, mut px2, mut px3] = base;
        let [mut mx0, mut mx1, mut mx2, mut mx3] = base;
        match j {
            0 => { px0 += eps; mx0 -= eps; }
            1 => { px1 += eps; mx1 -= eps; }
            2 => { px2 += eps; mx2 -= eps; }
            _ => { px3 += eps; mx3 -= eps; }
        }
        yoshida4_transverse_step(&mut px0, &mut px1, &mut px2, &mut px3, ds, k1, k2);
        yoshida4_transverse_step(&mut mx0, &mut mx1, &mut mx2, &mut mx3, ds, k1, k2);
        let plus = [px0, px1, px2, px3];
        let minus = [mx0, mx1, mx2, mx3];
        for i in 0..4 {
            jac[i][j] = (plus[i] - minus[i]) / (2.0 * eps);
        }
    }

    let det = det_4x4(jac);
    assert!((det - 1.0).abs() < 1.0e-10, "Symplectic det = {det}");
}

#[test]
fn sprint3_fodo_stable_vs_unstable_lattice() {
    // 4-meter FODO cell: Quad (0.2m), Drift (1.8m), Quad (0.2m), Drift (1.8m)
    // Stability condition: f > L/4 = 1.0 m  =>  k1 * l_q < 1.0  =>  k1 < 5.0 m^-2
    let build_fodo = |k1: f64, rigidity: f64| {
        let g = k1 * rigidity;
        Lattice::new(
            vec![
                LatticeElement::Quadrupole { length: 0.2, gradient: g },
                LatticeElement::Drift { length: 1.8 },
                LatticeElement::Quadrupole { length: 0.2, gradient: -g },
                LatticeElement::Drift { length: 1.8 },
            ],
            5.0,
        )
    };

    let rigidity = 1.0;
    let stable_lattice = build_fodo(2.5, rigidity); // f ~ 2.0 m > 1.0 m
    let (mut x, mut px, mut y, mut py) = (1.0e-3, 0.0, 1.0e-3, 0.0);
    let mut max_amp = 0.0_f64;
    for _ in 0..1_000 {
        stable_lattice.track_yoshida4_pass(&mut x, &mut px, &mut y, &mut py, rigidity, 8);
        max_amp = max_amp.max(x.abs()).max(y.abs());
    }
    assert!(max_amp < 0.01, "Stable FODO blew up: max_amp = {max_amp}");

    let unstable_lattice = build_fodo(6.5, rigidity); // f ~ 0.77 m < 1.0 m
    let (mut ux, mut upx, mut uy, mut upy) = (1.0e-3, 0.0, 1.0e-3, 0.0);
    for _ in 0..40 {
        unstable_lattice.track_yoshida4_pass(&mut ux, &mut upx, &mut uy, &mut upy, rigidity, 8);
    }
    assert!(ux.abs() > 1.0, "Unstable FODO failed to diverge: ux = {ux}");
}

#[test]
fn sprint4_rf_cavity_energy_gain() {
    let mut p = Particle::from_kinetic_energy_mev(
        ParticleSpecies::Proton,
        Vec3::zero(),
        Vec3::new(0.0, 1.0, 0.0),
        50.0,
    );
    let ek_before = p.kinetic_energy_ev();
    let v_rf = 500_000.0; // 500 kV
    let gap_len = 0.20;
    let phase = PI / 6.0; // sin(pi/6) = 0.5 => expected gain 250 keV
    let field = EmField::electric(Vec3::new(0.0, (v_rf / gap_len) * phase.sin(), 0.0));

    let v = p.velocity().norm();
    let dt = (gap_len / v) / 1_000.0;
    while p.pos.y < gap_len {
        let rem = gap_len - p.pos.y;
        let step_dt = dt.min(rem / p.velocity().y.max(1.0));
        boris_step(&mut p, field, step_dt);
        if (p.pos.y - gap_len).abs() < 1.0e-12 {
            break;
        }
    }
    let gain_ev = p.kinetic_energy_ev() - ek_before;
    let expected_ev = v_rf * phase.sin();
    let rel_err = (gain_ev - expected_ev).abs() / expected_ev;
    assert!(rel_err < 1.0e-3, "RF energy gain err = {rel_err}, gain = {gain_ev}");
}

#[test]
fn sprint4_isochronous_vs_classical_cyclotron() {
    let mut iso = AcceleratorSimulator::new(
        SimMode::IsochronousCyclotron,
        ParticleSpecies::Proton,
        8,
        2.0,
        1.5,
        0.0,
        1_200.0,
    );
    let mut cla = AcceleratorSimulator::new(
        SimMode::ClassicalCyclotron,
        ParticleSpecies::Proton,
        8,
        2.0,
        1.5,
        0.0,
        1_200.0,
    );
    iso.step(5_000);
    cla.step(5_000);

    assert!(
        iso.diagnostics.mean_ek_mev > 25.0,
        "Isochronous Ek too low: {}",
        iso.diagnostics.mean_ek_mev
    );
    assert!(
        iso.diagnostics.mean_ek_mev > cla.diagnostics.mean_ek_mev + 5.0,
        "Expected Isochronous ({:.2} MeV) to outperform Classical ({:.2} MeV) due to relativistic phase slip",
        iso.diagnostics.mean_ek_mev,
        cla.diagnostics.mean_ek_mev
    );
}

#[test]
fn sprint5_twiss_identity_and_emittance() {
    let bunch = BeamBunch::generate_ring_bunch(
        ParticleSpecies::Proton,
        128,
        5.0,
        100.0,
        0.005,
        0.001,
        0.003,
        0.0008,
        0.002,
        12345,
    );
    let diag = BeamDiagnostics::compute(&bunch, 5.0, 1.5);
    assert!(diag.twiss_x.emittance_geom > 0.0);
    assert!(diag.twiss_x.emittance_norm > diag.twiss_x.emittance_geom);
    assert!(
        diag.twiss_x.determinant_residual() < 1.0e-12,
        "Twiss identity residual = {}",
        diag.twiss_x.determinant_residual()
    );
}

#[test]
fn sprint5_liouville_emittance_conservation_fodo() {
    let fodo = Lattice::new(
        vec![
            LatticeElement::Quadrupole { length: 0.25, gradient: 2.0 },
            LatticeElement::Drift { length: 1.75 },
            LatticeElement::Quadrupole { length: 0.25, gradient: -2.0 },
            LatticeElement::Drift { length: 1.75 },
        ],
        5.0,
    );
    let mut xs = Vec::with_capacity(64);
    let mut pxs = Vec::with_capacity(64);
    let mut ys = Vec::with_capacity(64);
    let mut pys = Vec::with_capacity(64);

    let bunch = BeamBunch::generate_ring_bunch(
        ParticleSpecies::Proton,
        64,
        5.0,
        200.0,
        0.004,
        0.001,
        0.003,
        0.001,
        0.0,
        999,
    );
    for p in &bunch.particles {
        let r = (p.pos.x * p.pos.x + p.pos.y * p.pos.y).sqrt();
        xs.push(r - 5.0);
        pxs.push(p.u.x / p.u.norm());
        ys.push(p.pos.z);
        pys.push(p.u.z / p.u.norm());
    }
    let twiss_0 = TwissParameters::from_samples(&xs, &pxs, 1.0);

    for _ in 0..200 {
        for i in 0..64 {
            fodo.track_yoshida4_pass(&mut xs[i], &mut pxs[i], &mut ys[i], &mut pys[i], 1.5, 4);
        }
    }
    let twiss_final = TwissParameters::from_samples(&xs, &pxs, 1.0);
    let rel_eps_err =
        (twiss_final.emittance_geom - twiss_0.emittance_geom).abs() / twiss_0.emittance_geom;
    assert!(
        rel_eps_err < 1.0e-10,
        "Liouville emittance drift = {rel_eps_err}"
    );
}

#[test]
fn sprint5_wasm_ffi_buffers_valid() {
    sim_init(0, 0, 64, 150.0, 1.5, 2.5, 200.0);
    sim_step(20);
    assert_eq!(sim_get_particle_count(), 64);
    assert!(!sim_get_particles_ptr().is_null());
    assert!(!sim_get_phase_space_ptr().is_null());
    assert!(!sim_get_diagnostics_ptr().is_null());

    let mut sim = AcceleratorSimulator::new(
        SimMode::SynchrotronFodo,
        ParticleSpecies::Proton,
        64,
        150.0,
        1.5,
        2.5,
        200.0,
    );
    sim.step(20);
    assert!(sim.diagnostics_buf[0] > 1.0);
    assert!(sim.diagnostics_buf[2] > 140.0);
    assert!(sim.diagnostics_buf[18] < 1.0e-12);
}
