use eframe::egui;

use crate::simulation::affine_particle_in_cell_solver::AffineParticleInCellSolver;

pub struct FluidSolverApplication {
    affine_particle_in_cell_solver: AffineParticleInCellSolver,
    simulation_steps_per_frame: usize,
    time_step_seconds: f32,
}

impl FluidSolverApplication {
    pub fn new(_creation_context: &eframe::CreationContext<'_>) -> Self {
        Self {
            affine_particle_in_cell_solver: AffineParticleInCellSolver::new(48, 32, 0.03),
            simulation_steps_per_frame: 2,
            time_step_seconds: 1.0 / 120.0,
        }
    }
}

impl eframe::App for FluidSolverApplication {
    fn update(&mut self, context: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("controls_panel").show(context, |user_interface| {
            user_interface.horizontal(|user_interface| {
                user_interface.label("Simulation steps per frame");
                user_interface.add(
                    egui::Slider::new(&mut self.simulation_steps_per_frame, 1..=12).text("steps"),
                );
                user_interface.label("Time step (seconds)");
                user_interface.add(
                    egui::Slider::new(&mut self.time_step_seconds, 1.0 / 300.0..=1.0 / 30.0)
                        .logarithmic(true),
                );
            });
        });

        for _ in 0..self.simulation_steps_per_frame {
            self.affine_particle_in_cell_solver
                .advance(self.time_step_seconds);
        }

        egui::CentralPanel::default().show(context, |user_interface| {
            let available_space = user_interface.available_size();
            let (response, painter) =
                user_interface.allocate_painter(available_space, egui::Sense::hover());

            let simulation_size = self
                .affine_particle_in_cell_solver
                .simulation_dimensions_meters();

            let scaling_x = response.rect.width() / simulation_size.x;
            let scaling_y = response.rect.height() / simulation_size.y;
            let rendering_scale = scaling_x.min(scaling_y);

            let rendered_width = simulation_size.x * rendering_scale;
            let rendered_height = simulation_size.y * rendering_scale;
            let drawing_origin = egui::pos2(
                response.rect.left() + (response.rect.width() - rendered_width) * 0.5,
                response.rect.top() + (response.rect.height() - rendered_height) * 0.5,
            );

            painter.rect_stroke(
                egui::Rect::from_min_size(
                    drawing_origin,
                    egui::vec2(rendered_width, rendered_height),
                ),
                0.0,
                egui::Stroke::new(1.0, egui::Color32::LIGHT_GRAY),
            );

            for particle in self.affine_particle_in_cell_solver.particles() {
                let particle_screen_position = egui::pos2(
                    drawing_origin.x + particle.position.x * rendering_scale,
                    drawing_origin.y + particle.position.y * rendering_scale,
                );
                painter.circle_filled(
                    particle_screen_position,
                    2.0,
                    egui::Color32::from_rgb(30, 144, 255),
                );
            }
        });

        context.request_repaint();
    }
}
