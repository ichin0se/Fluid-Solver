use crate::simulation::grid::Grid;
use crate::simulation::particle::Particle;
use crate::simulation::vector2::Vector2;

pub struct AffineParticleInCellSolver {
    grid: Grid,
    particles: Vec<Particle>,
    gravity_meters_per_second_squared: Vector2,
}

impl AffineParticleInCellSolver {
    pub fn new(cell_count_x: usize, cell_count_y: usize, cell_size_meters: f32) -> Self {
        let node_count_x = cell_count_x + 1;
        let node_count_y = cell_count_y + 1;
        let mut particles = Vec::new();

        for particle_index_y in 0..(cell_count_y / 2) {
            for particle_index_x in 0..(cell_count_x / 2) {
                let position = Vector2::new(
                    (particle_index_x as f32 + 0.75) * cell_size_meters,
                    (particle_index_y as f32 + 0.75) * cell_size_meters,
                );
                particles.push(Particle::new(position));
            }
        }

        Self {
            grid: Grid::new(node_count_x, node_count_y, cell_size_meters),
            particles,
            gravity_meters_per_second_squared: Vector2::new(0.0, 9.81),
        }
    }

    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }

    pub fn simulation_dimensions_meters(&self) -> Vector2 {
        self.grid.dimensions_meters()
    }

    pub fn advance(&mut self, time_step_seconds: f32) {
        self.transfer_particle_state_to_grid();
        self.apply_external_forces(time_step_seconds);
        self.enforce_boundary_conditions();
        self.project_velocity_field_to_incompressible(30);
        self.transfer_grid_state_to_particles();
        self.advect_particles(time_step_seconds);
    }

    fn transfer_particle_state_to_grid(&mut self) {
        self.grid.reset_accumulators();

        for particle in &self.particles {
            let base_node_x = (particle.position.x / self.grid.cell_size_meters).floor() as isize;
            let base_node_y = (particle.position.y / self.grid.cell_size_meters).floor() as isize;

            for offset_y in 0..=1 {
                for offset_x in 0..=1 {
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
                    let node_index = self.grid.node_index(node_x_usize, node_y_usize);
                    let node_position =
                        self.grid.world_position_of_node(node_x_usize, node_y_usize);
                    let relative_position = node_position - particle.position;

                    let interpolation_weight = bilinear_weight(
                        particle.position,
                        node_position,
                        self.grid.cell_size_meters,
                    );

                    if interpolation_weight <= 0.0 {
                        continue;
                    }

                    let affine_velocity_component = Vector2::new(
                        particle.affine_velocity_matrix[0][0] * relative_position.x
                            + particle.affine_velocity_matrix[0][1] * relative_position.y,
                        particle.affine_velocity_matrix[1][0] * relative_position.x
                            + particle.affine_velocity_matrix[1][1] * relative_position.y,
                    );

                    let momentum_velocity = particle.velocity + affine_velocity_component;
                    let weighted_mass = interpolation_weight * particle.mass_kilograms;
                    self.grid.node_mass[node_index] += weighted_mass;
                    self.grid.node_velocity[node_index] += momentum_velocity * weighted_mass;
                }
            }
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

    fn project_velocity_field_to_incompressible(&mut self, solver_iteration_count: usize) {
        self.compute_divergence();
        self.grid.node_pressure.fill(0.0);

        let inverse_cell_size_squared =
            1.0 / (self.grid.cell_size_meters * self.grid.cell_size_meters);

        for _ in 0..solver_iteration_count {
            let previous_pressure_field = self.grid.node_pressure.clone();
            for node_y in 1..(self.grid.node_count_y - 1) {
                for node_x in 1..(self.grid.node_count_x - 1) {
                    let node_index = self.grid.node_index(node_x, node_y);
                    if self.grid.node_mass[node_index] <= 0.0 {
                        continue;
                    }

                    let left_pressure =
                        previous_pressure_field[self.grid.node_index(node_x - 1, node_y)];
                    let right_pressure =
                        previous_pressure_field[self.grid.node_index(node_x + 1, node_y)];
                    let bottom_pressure =
                        previous_pressure_field[self.grid.node_index(node_x, node_y - 1)];
                    let top_pressure =
                        previous_pressure_field[self.grid.node_index(node_x, node_y + 1)];

                    let divergence = self.grid.node_divergence[node_index];
                    self.grid.node_pressure[node_index] =
                        (left_pressure + right_pressure + bottom_pressure + top_pressure
                            - divergence / inverse_cell_size_squared)
                            * 0.25;
                }
            }
        }

        for node_y in 1..(self.grid.node_count_y - 1) {
            for node_x in 1..(self.grid.node_count_x - 1) {
                let node_index = self.grid.node_index(node_x, node_y);
                if self.grid.node_mass[node_index] <= 0.0 {
                    continue;
                }

                let pressure_gradient_x = (self.grid.node_pressure
                    [self.grid.node_index(node_x + 1, node_y)]
                    - self.grid.node_pressure[self.grid.node_index(node_x - 1, node_y)])
                    / (2.0 * self.grid.cell_size_meters);
                let pressure_gradient_y = (self.grid.node_pressure
                    [self.grid.node_index(node_x, node_y + 1)]
                    - self.grid.node_pressure[self.grid.node_index(node_x, node_y - 1)])
                    / (2.0 * self.grid.cell_size_meters);

                self.grid.node_velocity[node_index] -=
                    Vector2::new(pressure_gradient_x, pressure_gradient_y);
            }
        }
    }

    fn compute_divergence(&mut self) {
        self.grid.node_divergence.fill(0.0);

        for node_y in 1..(self.grid.node_count_y - 1) {
            for node_x in 1..(self.grid.node_count_x - 1) {
                let node_index = self.grid.node_index(node_x, node_y);
                if self.grid.node_mass[node_index] <= 0.0 {
                    continue;
                }

                let right_velocity =
                    self.grid.node_velocity[self.grid.node_index(node_x + 1, node_y)];
                let left_velocity =
                    self.grid.node_velocity[self.grid.node_index(node_x - 1, node_y)];
                let top_velocity =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y + 1)];
                let bottom_velocity =
                    self.grid.node_velocity[self.grid.node_index(node_x, node_y - 1)];

                self.grid.node_divergence[node_index] =
                    (right_velocity.x - left_velocity.x + top_velocity.y - bottom_velocity.y)
                        / (2.0 * self.grid.cell_size_meters);
            }
        }
    }

    fn transfer_grid_state_to_particles(&mut self) {
        for particle in &mut self.particles {
            let base_node_x = (particle.position.x / self.grid.cell_size_meters).floor() as isize;
            let base_node_y = (particle.position.y / self.grid.cell_size_meters).floor() as isize;

            let mut interpolated_velocity = Vector2::zero();
            let mut affine_numerator = [[0.0; 2]; 2];

            for offset_y in 0..=1 {
                for offset_x in 0..=1 {
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
                    let node_index = self.grid.node_index(node_x_usize, node_y_usize);
                    let node_position =
                        self.grid.world_position_of_node(node_x_usize, node_y_usize);
                    let relative_position = node_position - particle.position;

                    let interpolation_weight = bilinear_weight(
                        particle.position,
                        node_position,
                        self.grid.cell_size_meters,
                    );
                    if interpolation_weight <= 0.0 {
                        continue;
                    }

                    let node_velocity = self.grid.node_velocity[node_index];
                    interpolated_velocity += node_velocity * interpolation_weight;

                    affine_numerator[0][0] +=
                        interpolation_weight * node_velocity.x * relative_position.x;
                    affine_numerator[0][1] +=
                        interpolation_weight * node_velocity.x * relative_position.y;
                    affine_numerator[1][0] +=
                        interpolation_weight * node_velocity.y * relative_position.x;
                    affine_numerator[1][1] +=
                        interpolation_weight * node_velocity.y * relative_position.y;
                }
            }

            let inverse_moment = 4.0 / (self.grid.cell_size_meters * self.grid.cell_size_meters);
            particle.velocity = interpolated_velocity;
            particle.affine_velocity_matrix = [
                [
                    affine_numerator[0][0] * inverse_moment,
                    affine_numerator[0][1] * inverse_moment,
                ],
                [
                    affine_numerator[1][0] * inverse_moment,
                    affine_numerator[1][1] * inverse_moment,
                ],
            ];
        }
    }

    fn advect_particles(&mut self, time_step_seconds: f32) {
        let simulation_dimensions_meters = self.grid.dimensions_meters();
        let minimum_position = self.grid.cell_size_meters * 0.5;

        for particle in &mut self.particles {
            particle.position += particle.velocity * time_step_seconds;

            if particle.position.x < minimum_position {
                particle.position.x = minimum_position;
                particle.velocity.x = 0.0;
            }
            if particle.position.y < minimum_position {
                particle.position.y = minimum_position;
                particle.velocity.y = 0.0;
            }
            if particle.position.x > simulation_dimensions_meters.x - minimum_position {
                particle.position.x = simulation_dimensions_meters.x - minimum_position;
                particle.velocity.x = 0.0;
            }
            if particle.position.y > simulation_dimensions_meters.y - minimum_position {
                particle.position.y = simulation_dimensions_meters.y - minimum_position;
                particle.velocity.y = 0.0;
            }
        }
    }
}

fn bilinear_weight(
    sample_position: Vector2,
    lattice_node_position: Vector2,
    cell_size_meters: f32,
) -> f32 {
    let absolute_distance_x =
        ((sample_position.x - lattice_node_position.x) / cell_size_meters).abs();
    let absolute_distance_y =
        ((sample_position.y - lattice_node_position.y) / cell_size_meters).abs();

    if absolute_distance_x >= 1.0 || absolute_distance_y >= 1.0 {
        return 0.0;
    }

    (1.0 - absolute_distance_x) * (1.0 - absolute_distance_y)
}

#[cfg(test)]
mod tests {
    use super::{AffineParticleInCellSolver, bilinear_weight};
    use crate::simulation::vector2::Vector2;

    #[test]
    fn bilinear_weights_sum_to_one_inside_a_cell() {
        let sample_position = Vector2::new(0.25, 0.75);
        let cell_size_meters = 1.0;
        let nodes = [
            Vector2::new(0.0, 0.0),
            Vector2::new(1.0, 0.0),
            Vector2::new(0.0, 1.0),
            Vector2::new(1.0, 1.0),
        ];

        let weight_sum: f32 = nodes
            .iter()
            .map(|node| bilinear_weight(sample_position, *node, cell_size_meters))
            .sum();

        assert!((weight_sum - 1.0).abs() < 1e-6);
    }

    #[test]
    fn particle_positions_remain_inside_simulation_domain_after_advancing() {
        let mut solver = AffineParticleInCellSolver::new(16, 12, 0.05);
        let simulation_dimensions = solver.simulation_dimensions_meters();

        for _ in 0..120 {
            solver.advance(1.0 / 120.0);
        }

        for particle in solver.particles() {
            assert!(particle.position.x >= 0.0);
            assert!(particle.position.y >= 0.0);
            assert!(particle.position.x <= simulation_dimensions.x);
            assert!(particle.position.y <= simulation_dimensions.y);
        }
    }
}
