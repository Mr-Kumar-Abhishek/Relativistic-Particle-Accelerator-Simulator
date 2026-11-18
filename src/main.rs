//! Native CLI validation and physics benchmark runner for `particle-accelarator`.

use particle_accelarator::fields::EmField;
use particle_accelarator::integrator::boris_step;
use particle_accelarator::particle::{Particle, ParticleSpecies};
use particle_accelarator::simulator::{AcceleratorSimulator, SimMode};
use particle_accelarator::vec3::Vec3;

fn main() {
    println!("==================================================================");
    println!(" Relativistic Particle Accelerator Simulator (f64 Double Precision)");
    println!("==================================================================");

    // 1. 10,000-step Boris Energy Conservation Benchmark
    let mut probe = Particle::from_kinetic_energy_mev(
        ParticleSpecies::Proton,
        Vec3::zero(),
        Vec3::new(1.0, 0.0, 0.0),
        100.0,
    );
    let gamma0 = probe.gamma();
    let mag_field = EmField::magnetic(Vec3::new(0.0, 0.0, 2.0));
    for _ in 0..10_000 {
        boris_step(&mut probe, mag_field, 1.0e-11);
    }
    let rel_drift = (probe.gamma() - gamma0).abs() / gamma0;
    println!("\n[Benchmark 1] Relativistic Boris Pusher (10,000 steps in 2.0 T B-field)");
    println!("  Initial gamma       : {:.12}", gamma0);
    println!("  Final gamma         : {:.12}", probe.gamma());
    println!("  Relative drift      : {:.3e}  (target < 1.0e-14)", rel_drift);

    // 2. Synchrotron FODO Storage Ring Simulation
    let mut syn = AcceleratorSimulator::new(
        SimMode::SynchrotronFodo,
        ParticleSpecies::Proton,
        128,
        150.0,
        1.5,
        2.5,
        0.0, // Pure coasting ring first to verify emittance and energy conservation
    );
    let eps0 = syn.diagnostics.twiss_x.emittance_geom;
    syn.step(1_200); // 5 turns
    let d_syn = &syn.diagnostics;
    println!("\n[Benchmark 2] Synchrotron FODO Storage Ring (128 Protons @ 150 MeV, 5 turns)");
    println!("  Mean Ek / gamma / B : {:.4} MeV | {:.6} | {:.6}", d_syn.mean_ek_mev, d_syn.mean_gamma, d_syn.mean_beta);
    println!("  Magnetic Rigidity   : {:.6} T*m", d_syn.mean_rigidity_tm);
    println!("  RMS Envelope (x, y) : ({:.4} mm, {:.4} mm)", d_syn.twiss_x.rms_pos * 1.0e3, d_syn.twiss_y.rms_pos * 1.0e3);
    println!("  Geom Emittance x    : initial {:.4e} -> final {:.4e} m*rad", eps0, d_syn.twiss_x.emittance_geom);
    println!("  Twiss (alpha, beta) : ({:.4}, {:.4} m) | |bg - a^2 - 1| = {:.3e}", d_syn.twiss_x.alpha, d_syn.twiss_x.beta, d_syn.twiss_x.determinant_residual());
    println!("  Energy Conservation : {:.3e}", d_syn.energy_drift_rel);

    // 3. Isochronous vs. Classical Cyclotron Acceleration Comparison
    let mut iso = AcceleratorSimulator::new(
        SimMode::IsochronousCyclotron,
        ParticleSpecies::Proton,
        16,
        2.0,
        1.5,
        0.0,
        1_200.0,
    );
    let mut cla = AcceleratorSimulator::new(
        SimMode::ClassicalCyclotron,
        ParticleSpecies::Proton,
        16,
        2.0,
        1.5,
        0.0,
        1_200.0,
    );
    iso.step(5_000);
    cla.step(5_000);
    println!("\n[Benchmark 3] Cyclotron RF Acceleration Comparison (5,000 steps, 1.2 MV Dee)");
    println!(
        "  Isochronous Cyclotron : Ek = {:.3} MeV (gamma = {:.5}, turns = {:.1})",
        iso.diagnostics.mean_ek_mev, iso.diagnostics.mean_gamma, iso.turns_completed
    );
    println!(
        "  Classical Cyclotron   : Ek = {:.3} MeV (gamma = {:.5}, turns = {:.1})",
        cla.diagnostics.mean_ek_mev, cla.diagnostics.mean_gamma, cla.turns_completed
    );
    println!("\nAll CLI validation benchmarks completed successfully.");
}
