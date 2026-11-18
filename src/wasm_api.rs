//! Zero-copy C-ABI WebAssembly interface for real-time browser execution.

use core::cell::UnsafeCell;

use crate::particle::ParticleSpecies;
use crate::simulator::{AcceleratorSimulator, SimMode};

struct GlobalSim(UnsafeCell<Option<AcceleratorSimulator>>);

// WebAssembly executes single-threaded unless SharedArrayBuffer threads are explicitly enabled.
unsafe impl Sync for GlobalSim {}

static SIM_INSTANCE: GlobalSim = GlobalSim(UnsafeCell::new(None));

fn with_sim<R>(f: impl FnOnce(&mut AcceleratorSimulator) -> R, default: R) -> R {
    unsafe {
        if let Some(sim) = (*SIM_INSTANCE.0.get()).as_mut() {
            f(sim)
        } else {
            default
        }
    }
}

/// Initializes (or re-initializes) the global simulator instance inside WASM linear memory.
#[cfg_attr(target_arch = "wasm32", no_mangle)]
pub extern "C" fn sim_init(
    mode: u32,
    species_id: u32,
    particle_count: u32,
    energy_mev: f64,
    b_field_t: f64,
    quad_gradient_tm: f64,
    rf_voltage_kv: f64,
) {
    let sim = AcceleratorSimulator::new(
        SimMode::from_id(mode),
        ParticleSpecies::from_id(species_id),
        particle_count as usize,
        energy_mev,
        b_field_t,
        quad_gradient_tm,
        rf_voltage_kv,
    );
    unsafe {
        *SIM_INSTANCE.0.get() = Some(sim);
    }
}

/// Advances the active simulation by `substeps` time steps.
#[cfg_attr(target_arch = "wasm32", no_mangle)]
pub extern "C" fn sim_step(substeps: u32) {
    with_sim(|sim| sim.step(substeps as usize), ());
}

/// Returns the linear-memory pointer to the flat `[x, y, z, vx, vy, vz, gamma, ek_mev]` particle array.
#[cfg_attr(target_arch = "wasm32", no_mangle)]
pub extern "C" fn sim_get_particles_ptr() -> *const f64 {
    with_sim(|sim| sim.particles_buf.as_ptr(), core::ptr::null())
}

/// Returns the linear-memory pointer to the flat `[x, xp, y, yp, dphi, de_mev]` phase-space array.
#[cfg_attr(target_arch = "wasm32", no_mangle)]
pub extern "C" fn sim_get_phase_space_ptr() -> *const f64 {
    with_sim(|sim| sim.phase_space_buf.as_ptr(), core::ptr::null())
}

/// Returns the linear-memory pointer to the 20-element diagnostics telemetry array.
#[cfg_attr(target_arch = "wasm32", no_mangle)]
pub extern "C" fn sim_get_diagnostics_ptr() -> *const f64 {
    with_sim(|sim| sim.diagnostics_buf.as_ptr(), core::ptr::null())
}

/// Returns the number of particles currently allocated in the bunch.
#[cfg_attr(target_arch = "wasm32", no_mangle)]
pub extern "C" fn sim_get_particle_count() -> u32 {
    with_sim(|sim| sim.bunch.particles.len() as u32, 0)
}

/// Returns the design ring or cyclotron extraction radius (m).
#[cfg_attr(target_arch = "wasm32", no_mangle)]
pub extern "C" fn sim_get_ring_radius() -> f64 {
    with_sim(|sim| sim.diagnostics_buf[19], 1.0)
}
