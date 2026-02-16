use crate::simulation::grid::Grid;
use crate::simulation::particle::Particle;
use crate::simulation::vector2::Vector2;

/// APIC/FLIP solver inspired by commonly used production pipelines:
/// - Jiang et al., "The Affine Particle-In-Cell Method" (SIGGRAPH 2015).
/// - Zhu and Bridson, "Animating Sand as a Fluid" (SIGGRAPH 2005, FLIP variant usage).
/// - Bridson, "Fluid Simulation for Computer Graphics" (pressure projection practice).
pub struct AffineParticleInCellSolver {
    grid: Grid,
    previous_grid_velocity: Vec<Vector2>,
    particles: Vec<Particle>,
    gravity_meters_per_second_squared: Vector2,
    fluid_density_kilograms_per_cubic_meter: f32,
    pressure_solver_iteration_count: usize,
    velocity_pic_blend: f32,
    collision_restitution: f32,
    collision_tangential_damping: f32,
    particle_collision_radius_meters: f32,
    vorticity_confinement_strength: f32,
}

impl AffineParticleInCellSolver {
    pub fn new(cell_count_x: usize, cell_count_y: usize, cell_size_meters: f32) -> Self {
        let node_count_x = cell_count_x + 1;
        let node_count_y = cell_count_y + 1;
        let grid = Grid::new(node_count_x, node_count_y, cell_size_meters);

        let mut solver = Self {
            previous_grid_velocity: vec![Vector2::zero(); node_count_x * node_count_y],
            grid,
            particles: Vec::new(),
            gravity_meters_per_second_squared: Vector2::new(0.0, 9.81),
            fluid_density_kilograms_per_cubic_meter: 1000.0,
            pressure_solver_iteration_count: 120,
            velocity_pic_blend: 0.05,
            collision_restitution: 0.45,
            collision_tangential_damping: 0.02,
            particle_collision_radius_meters: cell_size_meters * 0.35,
            vorticity_confinement_strength: 0.45,
        };
        solver.reset();
        solver
    }

    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }

    pub fn simulation_dimensions_meters(&self) -> Vector2 {
        self.grid.dimensions_meters()
    }

    pub fn velocity_pic_blend(&self) -> f32 {
        self.velocity_pic_blend
    }

    pub fn set_velocity_pic_blend(&mut self, velocity_pic_blend: f32) {
        self.velocity_pic_blend = velocity_pic_blend.clamp(0.0, 1.0);
    }

    pub fn reset(&mut self) {
        self.particles = self.create_initial_particles();
        self.grid.node_velocity.fill(Vector2::zero());
        self.grid.node_mass.fill(0.0);
        self.grid.node_pressure.fill(0.0);
        self.grid.node_divergence.fill(0.0);
        self.previous_grid_velocity.fill(Vector2::zero());
    }

    pub fn advance(&mut self, time_step_seconds: f32) {
        self.transfer_particle_state_to_grid();
        self.previous_grid_velocity
            .clone_from(&self.grid.node_velocity);

        self.apply_external_forces(time_step_seconds);
        self.apply_vorticity_confinement(time_step_seconds);
        self.enforce_boundary_conditions();

        self.project_velocity_field_to_incompressible(time_step_seconds);
        self.enforce_boundary_conditions();

        self.transfer_grid_state_to_particles();
        self.advect_particles(time_step_seconds);
    }

    fn create_initial_particles(&self) -> Vec<Particle> {
        let cell_count_x = self.grid.node_count_x - 1;
        let cell_count_y = self.grid.node_count_y - 1;

        let mut particles = Vec::new();
        for particle_index_y in 0..(cell_count_y * 3 / 5) {
            for particle_index_x in 0..(cell_count_x * 4 / 7) {
                let position = Vector2::new(
                    (particle_index_x as f32 + 0.5) * self.grid.cell_size_meters,
                    (particle_index_y as f32 + 0.5) * self.grid.cell_size_meters,
                );
                let mut particle = Particle::new(position);

                let horizontal_ratio =
                    (particle_index_x as f32 / cell_count_x.max(1) as f32).clamp(0.0, 1.0);
                particle.velocity.x = (horizontal_ratio - 0.5) * 1.0;

                particles.push(particle);
            }
        }

        particles
    }

    fn transfer_particle_state_to_grid(&mut self) {
        self.grid.reset_accumulators();

        for particle in &self.particles {
            self.for_each_neighbor_node(particle.position, |node_index, node_position, weight| {
                let relative_position = node_position - particle.position;

                let affine_velocity_component = Vector2::new(
                    particle.affine_velocity_matrix[0][0] * relative_position.x
                        + particle.affine_velocity_matrix[0][1] * relative_position.y,
                    particle.affine_velocity_matrix[1][0] * relative_position.x
                        + particle.affine_velocity_matrix[1][1] * relative_position.y,
                );

                let particle_velocity_with_affine = particle.velocity + affine_velocity_component;
                let weighted_mass = weight * particle.mass_kilograms;

                self.grid.node_mass[node_index] += weighted_mass;
                self.grid.node_velocity[node_index] +=
                    particle_velocity_with_affine * weighted_mass;
            });
        }

        for node_index in 0..self.grid.node_mass.len() {
            let node_mass = self.grid.node_mass[node_index];
            if node_mass > 0.0 {
                self.grid.node_velocity[node_index] =
                    self.grid.node_velocity[node_index] / node_mass;
            }
        }
    }

    fn apply_external_forces(&mut self, time_step_seconds: f32) {
        for node_index in 0..self.grid.node_mass.len() {
            if self.grid.node_mass[node_index] > 0.0 {
                self.grid.node_velocity[node_index] +=
                    self.gravity_meters_per_second_squared * time_step_seconds;
            }
        }
    }

    fn apply_vorticity_confinement(&mut self, time_step_seconds: f32) {
        let node_count_total = self.grid.node_count_x * self.grid.node_count_y;
        let mut curl_magnitude = vec![0.0; node_count_total];

        for node_y in 1..(self.grid.node_count_y - 1) {
            for node_x in 1..(self.grid.node_count_x - 1) {
                let node_index = self.grid.node_index(node_x, node_y);
                if self.grid.node_mass[node_index] <= 0.0 {
                    continue;
                }

                let velocity_right =
                    self.grid.node_velocity[self.grid.node_index(node_x + 1, node_y)];
                let velocity_left =
                    self.grid.node_velocity[self.grid.node_index(node_x - 1, node_y)];
                let velocity_top =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y + 1)];
                let velocity_bottom =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y - 1)];

                let curl =
                    (velocity_right.y - velocity_left.y - (velocity_top.x - velocity_bottom.x))
                        / (2.0 * self.grid.cell_size_meters);
                curl_magnitude[node_index] = curl.abs();
            }
        }

        for node_y in 2..(self.grid.node_count_y - 2) {
            for node_x in 2..(self.grid.node_count_x - 2) {
                let node_index = self.grid.node_index(node_x, node_y);
                if self.grid.node_mass[node_index] <= 0.0 {
                    continue;
                }

                let gradient_x = (curl_magnitude[self.grid.node_index(node_x + 1, node_y)]
                    - curl_magnitude[self.grid.node_index(node_x - 1, node_y)])
                    / (2.0 * self.grid.cell_size_meters);
                let gradient_y = (curl_magnitude[self.grid.node_index(node_x, node_y + 1)]
                    - curl_magnitude[self.grid.node_index(node_x, node_y - 1)])
                    / (2.0 * self.grid.cell_size_meters);

                let gradient_length = (gradient_x * gradient_x + gradient_y * gradient_y).sqrt();
                if gradient_length < 1e-5 {
                    continue;
                }

                let normalized_x = gradient_x / gradient_length;
                let normalized_y = gradient_y / gradient_length;

                let velocity_right =
                    self.grid.node_velocity[self.grid.node_index(node_x + 1, node_y)];
                let velocity_left =
                    self.grid.node_velocity[self.grid.node_index(node_x - 1, node_y)];
                let velocity_top =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y + 1)];
                let velocity_bottom =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y - 1)];
                let curl =
                    (velocity_right.y - velocity_left.y - (velocity_top.x - velocity_bottom.x))
                        / (2.0 * self.grid.cell_size_meters);

                let confinement_force = Vector2::new(normalized_y * curl, -normalized_x * curl)
                    * self.vorticity_confinement_strength;
                self.grid.node_velocity[node_index] += confinement_force * time_step_seconds;
            }
        }
    }

    fn enforce_boundary_conditions(&mut self) {
        for node_y in 0..self.grid.node_count_y {
            for node_x in 0..self.grid.node_count_x {
                if node_x == 0
                    || node_y == 0
                    || node_x == self.grid.node_count_x - 1
                    || node_y == self.grid.node_count_y - 1
                {
                    let node_index = self.grid.node_index(node_x, node_y);
                    self.grid.node_velocity[node_index] = Vector2::zero();
                }
            }
        }
    }

    fn project_velocity_field_to_incompressible(&mut self, time_step_seconds: f32) {
        self.compute_divergence();
        self.grid.node_pressure.fill(0.0);

        let cell_size = self.grid.cell_size_meters;
        let cell_size_squared = cell_size * cell_size;

        for _ in 0..self.pressure_solver_iteration_count {
            let previous_pressure_field = self.grid.node_pressure.clone();

            for node_y in 1..(self.grid.node_count_y - 1) {
                for node_x in 1..(self.grid.node_count_x - 1) {
                    let node_index = self.grid.node_index(node_x, node_y);
                    if self.grid.node_mass[node_index] <= 0.0 {
                        continue;
                    }

                    let mut pressure_neighbor_sum = 0.0;
                    let mut active_neighbor_count = 0.0;

                    let neighbors = [
                        (node_x - 1, node_y),
                        (node_x + 1, node_y),
                        (node_x, node_y - 1),
                        (node_x, node_y + 1),
                    ];

                    for (neighbor_x, neighbor_y) in neighbors {
                        let neighbor_index = self.grid.node_index(neighbor_x, neighbor_y);
                        if self.grid.node_mass[neighbor_index] > 0.0 {
                            pressure_neighbor_sum += previous_pressure_field[neighbor_index];
                            active_neighbor_count += 1.0;
                        }
                    }

                    if active_neighbor_count <= 0.0 {
                        continue;
                    }

                    let right_hand_side = self.fluid_density_kilograms_per_cubic_meter
                        / time_step_seconds
                        * self.grid.node_divergence[node_index];

                    self.grid.node_pressure[node_index] = (pressure_neighbor_sum
                        - right_hand_side * cell_size_squared)
                        / active_neighbor_count;
                }
            }
        }

        let pressure_gradient_scale =
            time_step_seconds / self.fluid_density_kilograms_per_cubic_meter;

        for node_y in 1..(self.grid.node_count_y - 1) {
            for node_x in 1..(self.grid.node_count_x - 1) {
                let node_index = self.grid.node_index(node_x, node_y);
                if self.grid.node_mass[node_index] <= 0.0 {
                    continue;
                }

                let pressure_gradient_x = (self.grid.node_pressure
                    [self.grid.node_index(node_x + 1, node_y)]
                    - self.grid.node_pressure[self.grid.node_index(node_x - 1, node_y)])
                    / (2.0 * cell_size);
                let pressure_gradient_y = (self.grid.node_pressure
                    [self.grid.node_index(node_x, node_y + 1)]
                    - self.grid.node_pressure[self.grid.node_index(node_x, node_y - 1)])
                    / (2.0 * cell_size);

                self.grid.node_velocity[node_index] -=
                    Vector2::new(pressure_gradient_x, pressure_gradient_y)
                        * pressure_gradient_scale;
            }
        }
    }

    fn compute_divergence(&mut self) {
        self.grid.node_divergence.fill(0.0);
        let inverse_double_cell_size = 1.0 / (2.0 * self.grid.cell_size_meters);

        for node_y in 1..(self.grid.node_count_y - 1) {
            for node_x in 1..(self.grid.node_count_x - 1) {
                let node_index = self.grid.node_index(node_x, node_y);
                if self.grid.node_mass[node_index] <= 0.0 {
                    continue;
                }

                let velocity_right =
                    self.grid.node_velocity[self.grid.node_index(node_x + 1, node_y)];
                let velocity_left =
                    self.grid.node_velocity[self.grid.node_index(node_x - 1, node_y)];
                let velocity_top =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y + 1)];
                let velocity_bottom =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y - 1)];

                self.grid.node_divergence[node_index] =
                    (velocity_right.x - velocity_left.x + velocity_top.y - velocity_bottom.y)
                        * inverse_double_cell_size;
            }
        }
    }

    fn transfer_grid_state_to_particles(&mut self) {
        for particle in &mut self.particles {
            let mut interpolated_velocity_from_current_grid = Vector2::zero();
            let mut interpolated_velocity_change_from_grid = Vector2::zero();
            let mut affine_numerator = [[0.0; 2]; 2];

            self.for_each_neighbor_node(particle.position, |node_index, node_position, weight| {
                let current_node_velocity = self.grid.node_velocity[node_index];
                let previous_node_velocity = self.previous_grid_velocity[node_index];
                let node_velocity_change = current_node_velocity - previous_node_velocity;
                let relative_position = node_position - particle.position;

                interpolated_velocity_from_current_grid += current_node_velocity * weight;
                interpolated_velocity_change_from_grid += node_velocity_change * weight;

                affine_numerator[0][0] += weight * current_node_velocity.x * relative_position.x;
                affine_numerator[0][1] += weight * current_node_velocity.x * relative_position.y;
                affine_numerator[1][0] += weight * current_node_velocity.y * relative_position.x;
                affine_numerator[1][1] += weight * current_node_velocity.y * relative_position.y;
            });

            let velocity_from_flip = particle.velocity + interpolated_velocity_change_from_grid;
            let velocity_from_pic = interpolated_velocity_from_current_grid;
            particle.velocity = velocity_from_flip * (1.0 - self.velocity_pic_blend)
                + velocity_from_pic * self.velocity_pic_blend;

            let inverse_second_moment =
                4.0 / (self.grid.cell_size_meters * self.grid.cell_size_meters);
            particle.affine_velocity_matrix = [
                [
                    affine_numerator[0][0] * inverse_second_moment,
                    affine_numerator[0][1] * inverse_second_moment,
                ],
                [
                    affine_numerator[1][0] * inverse_second_moment,
                    affine_numerator[1][1] * inverse_second_moment,
                ],
            ];
        }
    }

    fn advect_particles(&mut self, time_step_seconds: f32) {
        let minimum_position = self.particle_collision_radius_meters;
        let simulation_dimensions_meters = self.grid.dimensions_meters();

        for particle in &mut self.particles {
            let position_stage_0 = particle.position;
            let velocity_stage_0 = self.sample_grid_velocity(position_stage_0, true);
            let position_stage_1 = position_stage_0 + velocity_stage_0 * time_step_seconds;

            let velocity_stage_1 = self.sample_grid_velocity(position_stage_1, true);
            let position_stage_2 = position_stage_0 * 0.75
                + (position_stage_1 + velocity_stage_1 * time_step_seconds) * 0.25;

            let velocity_stage_2 = self.sample_grid_velocity(position_stage_2, true);
            particle.position = position_stage_0 * (1.0 / 3.0)
                + (position_stage_2 + velocity_stage_2 * time_step_seconds) * (2.0 / 3.0);

            if particle.position.x < minimum_position {
                particle.position.x = minimum_position;
                particle.velocity.x = -particle.velocity.x * self.collision_restitution;
                particle.velocity.y *= 1.0 - self.collision_tangential_damping;
            }
            if particle.position.y < minimum_position {
                particle.position.y = minimum_position;
                particle.velocity.y = -particle.velocity.y * self.collision_restitution;
                particle.velocity.x *= 1.0 - self.collision_tangential_damping;
            }
            if particle.position.x > simulation_dimensions_meters.x - minimum_position {
                particle.position.x = simulation_dimensions_meters.x - minimum_position;
                particle.velocity.x = -particle.velocity.x * self.collision_restitution;
                particle.velocity.y *= 1.0 - self.collision_tangential_damping;
            }
            if particle.position.y > simulation_dimensions_meters.y - minimum_position {
                particle.position.y = simulation_dimensions_meters.y - minimum_position;
                particle.velocity.y = -particle.velocity.y * self.collision_restitution;
                particle.velocity.x *= 1.0 - self.collision_tangential_damping;
            }
        }
    }

    fn sample_grid_velocity(
        &self,
        sample_position: Vector2,
        use_current_velocity: bool,
    ) -> Vector2 {
        let mut sampled_velocity = Vector2::zero();

        self.for_each_neighbor_node(sample_position, |node_index, _node_position, weight| {
            let source_velocity = if use_current_velocity {
                self.grid.node_velocity[node_index]
            } else {
                self.previous_grid_velocity[node_index]
            };
            sampled_velocity += source_velocity * weight;
        });

        sampled_velocity
    }

    fn for_each_neighbor_node(
        &self,
        sample_position: Vector2,
        mut callback: impl FnMut(usize, Vector2, f32),
    ) {
        let base_node_x = (sample_position.x / self.grid.cell_size_meters).floor() as isize;
        let base_node_y = (sample_position.y / self.grid.cell_size_meters).floor() as isize;

        for offset_y in -1..=1 {
            for offset_x in -1..=1 {
                let node_x = base_node_x + offset_x;
                let node_y = base_node_y + offset_y;

                if node_x < 0
                    || node_y < 0
                    || node_x as usize >= self.grid.node_count_x
                    || node_y as usize >= self.grid.node_count_y
                {
                    continue;
                }

                let node_x_usize = node_x as usize;
                let node_y_usize = node_y as usize;
                let node_position = self.grid.world_position_of_node(node_x_usize, node_y_usize);

                let weight_x = quadratic_b_spline_weight(
                    (sample_position.x - node_position.x) / self.grid.cell_size_meters,
                );
                let weight_y = quadratic_b_spline_weight(
                    (sample_position.y - node_position.y) / self.grid.cell_size_meters,
                );
                let weight = weight_x * weight_y;

                if weight > 0.0 {
                    let node_index = self.grid.node_index(node_x_usize, node_y_usize);
                    callback(node_index, node_position, weight);
                }
            }
        }
    }
}

fn quadratic_b_spline_weight(normalized_distance: f32) -> f32 {
    let absolute_distance = normalized_distance.abs();

    if absolute_distance < 0.5 {
        0.75 - absolute_distance * absolute_distance
    } else if absolute_distance < 1.5 {
        let difference = 1.5 - absolute_distance;
        0.5 * difference * difference
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::{AffineParticleInCellSolver, quadratic_b_spline_weight};

    #[test]
    fn quadratic_b_spline_weights_are_positive_near_center() {
        assert!(quadratic_b_spline_weight(0.0) > 0.0);
        assert!(quadratic_b_spline_weight(0.8) > 0.0);
        assert_eq!(quadratic_b_spline_weight(2.0), 0.0);
    }

    #[test]
    fn particle_positions_remain_inside_simulation_domain_after_advancing() {
        let mut solver = AffineParticleInCellSolver::new(24, 18, 0.04);
        let simulation_dimensions = solver.simulation_dimensions_meters();

        for _ in 0..240 {
            solver.advance(1.0 / 120.0);
        }

        for particle in solver.particles() {
            assert!(particle.position.x >= 0.0);
            assert!(particle.position.y >= 0.0);
            assert!(particle.position.x <= simulation_dimensions.x);
            assert!(particle.position.y <= simulation_dimensions.y);
        }
    }

    #[test]
    fn reset_restores_initial_particle_count() {
        let mut solver = AffineParticleInCellSolver::new(20, 16, 0.04);
        let initial_particle_count = solver.particles().len();

        for _ in 0..30 {
            solver.advance(1.0 / 120.0);
        }
        solver.reset();

        assert_eq!(solver.particles().len(), initial_particle_count);
    }

    #[test]
    fn velocity_pic_blend_is_clamped() {
        let mut solver = AffineParticleInCellSolver::new(16, 12, 0.04);

        solver.set_velocity_pic_blend(-2.0);
        assert_eq!(solver.velocity_pic_blend(), 0.0);

        solver.set_velocity_pic_blend(2.0);
        assert_eq!(solver.velocity_pic_blend(), 1.0);
    }
}
