mod application;
mod simulation;

use application::fluid_solver_application::FluidSolverApplication;

fn main() -> eframe::Result<()> {
    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "APIC Liquid Solver",
        native_options,
        Box::new(|creation_context| Ok(Box::new(FluidSolverApplication::new(creation_context)))),
    )
}
